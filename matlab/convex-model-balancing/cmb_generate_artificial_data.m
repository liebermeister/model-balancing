function [kinetics, prior, bounds, data, true, kinetic_data, state_data] = cmb_generate_artificial_data(network, cmb_options, q_info, c_init)

% [kinetics, prior, bounds, data, true, kinetic_data, state_data] = cmb_generate_artificial_data(network, cmb_options, q_info)
%
% Generate an artificial data set for a given metabolic network structure
% 
%   network: network structure
%   cmb_options  struct containing options
%   cmb_options.ns:      number of steady states
%   q_info       struct describing the dependencies between model variables
%   c_init       (optional) matrix; initial guess of (non-logarithmic) metabolite concentrations; 
%                           also used for generating the artificial data 
%
%   prior        struct describing priors in the optimality problem (details see cmb_make_prior)
%   bounds       struct describing bounds in the optimality problem (details see cmb_make_bounds)
%   data         struct describing data used in the optimality problem
%   true         struct describing true model variables (optional; only for models with artificial data)

eval(default('c_init','[]'));

% --------------------------------------------------------------
%% Set global variables to speed up function modular_velocities
global global_structure_matrices 
global_structure_matrices = 1;
global Mplus Mminus Wplus Wminus nm nr ind_M ind_Wp ind_Wm
N = network.N; W = network.regulation_matrix; ind_ext = find(network.external); h = ones(size(N,2),1);
[Mplus, Mminus, Wplus, Wminus, nm, nr, N_int,ind_M,ind_Wp,ind_Wm] = make_structure_matrices(N,W,ind_ext,h);
% --------------------------------------------------------------

[nm,nr] = size(network.N);

if ~strcmp(cmb_options.ecm_score,'emc4cm'), error('wrong rate law'); end 

switch cmb_options.ecm_score
  case 'emc4cm',
    network.kinetics = set_kinetics(network,'cs');
  otherwise
    error('unsupported rate laws');
end

% --------------------------------------------------------------
% kinetic constants

rng(cmb_options.random_seed,'twister');

%% for simplicity: randomly choose std chem pot
mu0                  = cmb_options.quantities.mu0.std * randn(nm,1);
network.kinetics.Keq = exp(- network.N' * mu0/RT); 
network.kinetics.KV  = exp(gaussian_random_matrix(cmb_options.quantities.KV.mean_ln,cmb_options.quantities.KV.std_ln,log(1.1*cmb_options.quantities.KV.min),log(0.9*cmb_options.quantities.KV.max),nr,1));
network.kinetics.KM  = double(network.kinetics.KM~=0) .* exp(gaussian_random_matrix(cmb_options.quantities.KM.mean_ln,cmb_options.quantities.KM.std_ln,log(1.1*cmb_options.quantities.KM.min),log(0.9*cmb_options.quantities.KM.max),nr,nm));
if isempty(c_init),
  c_init = ones(length(network.metabolites),cmb_options.ns);
end
network.kinetics.c   = c_init;
network.kinetics.u   = network.kinetics.u .* exp(0.5*randn(nr,1));


% -------------------------------------------------------------

% [qind, qall, M_qind_to_qall, M_qall_to_qind] = cmb_analyse_network(network);
% NOTE THAT THE FUNCTION NEEDS TO BE FIXED
 
% -------------------------------------------------------------
% Initialise other variables

[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);

keq = network.kinetics.Keq;


% --------------------------------------------------------------
% Bounds

% FIX: rather get these numbers from parameter prior table

bounds = cmb_make_bounds(network,q_info,cmb_options);

% --------------------------------------------------------------
% Priors (kinetic constants)

prior = cmb_make_prior(network,q_info,cmb_options);


% -------------------------------------------------------------
% Generate "true" data points based on existing data set;
% use feasible state, randomise enzyme and ext. metab levels (multiplicative random variations),
% and find new steady states

clear true

ns = cmb_options.ns;

true.q = cmb_kinetics_to_q(network, cmb_options, q_info);

% --------------------------------------------------------------
% metabolic states

rng(cmb_options.random_seed,'twister');

% compute one feasible steady state ("reference" )
% factors 0.1 are used to obtain concentrations in reasonable ranges

e_ref = cmb_options.metabolic_prior_e_geom_mean * [cmb_options.metabolic_artificial_e_geom_std .^ randn(nr,1)];
network.kinetics.u = e_ref;

[c_ref, v_ref] = network_steady_state(network, network.kinetics.c);
x_ref  = log(c_ref);

% check correctness and driving forces
% pp   = cmb_make_pp(network);
% pp.v = v_ref;
% [~, e_ref] = ecm_get_score(cmb_options.ecm_score,x_ref,pp);
% [network.kinetics.u, e_ref]
% theta = driving_force(c_ref,v_ref,network)
% pm(RT * theta, network.actions)

% perturb the reference state to obtain several steady states ("samples")

X_init = repmat(x_ref,1,ns)  + log(cmb_options.metabolic_artificial_c_geom_std) * randn(nm,ns);
true.E = repmat(e_ref,1,ns) .*    [cmb_options.metabolic_artificial_e_geom_std .^ randn(nr,ns)];

for j=1:ns,
  my_network = network;
  my_network.kinetics.c = exp(X_init(:,j));
  my_network.kinetics.u = true.E(:,j);
  [c_ss,v_ss] = network_steady_state(my_network,my_network.kinetics.c);
  true.X(:,j) = log(c_ss);
  true.V(:,j) = v_ss;
  true.A_forward(:,j) = sign(v_ss) .* [log(network.kinetics.Keq) - network.N' * true.X(:,j)];
end

true.kinetics = network.kinetics;
[true.kinetics.Kcatf, true.kinetics.Kcatr] = modular_KV_Keq_to_kcat(network.N,network.kinetics);


% --------------------------------------------------------------
% Data matrices

for it = 1:ns,
  data.samples{it,1} = ['S' num2str(it)];
end

data.V.mean = true.V;
data.X.mean = true.X;
data.E.mean = true.E;

data.X.std  = log(cmb_options.data_C_geom_std)    * ones(size(data.X.mean));
data.E.std  =    [cmb_options.data_E_geom_std-1]  * data.E.mean      + 0.01 * max(max(abs(data.E.mean)));
data.V.std  =    [cmb_options.data_V_geom_std-1]  * abs(data.V.mean) + 0.01 * max(max(abs(data.V.mean)));

% noise in artificial data (NOTE: the resulting "states" may be thermo-physiologically infeasible!)

if cmb_options.use_artificial_noise,
  if cmb_options.verbose,
    display('Using artificial metabolic data with noise');
  end
  data.V.mean = data.V.mean + data.V.std .* randn(nr,ns);
  ind_wrong_sign = find(data.V.mean .* data.V.mean <0);
  data.V.mean(ind_wrong_sign) = -data.V.mean(ind_wrong_sign);
  data.X.mean = data.X.mean + data.X.std .* randn(nm,ns);
  data.E.mean = data.E.mean + data.E.std .* randn(nr,ns);
else
  if cmb_options.verbose,
    display('Using artificial metabolic data without noise');
  end
end

data.qall.mean = cmb_q_to_qall(true.q,q_info);
data.qall.std  = log(cmb_options.data_kin_geom_std) * ones(size(data.qall.mean));

if cmb_options.use_kinetic_data_noise,
  if cmb_options.verbose,
    display('Generating artificial kinetic data with noise');
  end
  data.qall.mean = data.qall.mean + data.qall.std .* randn(size(data.qall.mean));
else
  if cmb_options.verbose,
    display('Generating artificial kinetic data without noise');
  end
end


% --------------------------------------------------------------
% check whether true values satisfy bounds

flag_ok = 1;

if find(true.q < bounds.q_min),                              warning('Lower bounds for q violated'); flag_ok = 0; end
if find(true.q > bounds.q_max),                              warning('Upper bounds for q violated'); flag_ok = 0; end
if find(true.X < repmat(bounds.x_min,1,ns)),                 warning('Lower bounds for X violated'); flag_ok = 0; end
if find(true.X > repmat(bounds.x_max,1,ns)),                 warning('Upper bounds for X violated'); flag_ok = 0; end
if find(true.E < repmat(bounds.e_min,1,ns)),                 warning('Lower bounds for E violated'); flag_ok = 0; end
if find(true.E > repmat(bounds.e_max,1,ns)),                 warning('Upper bounds for E violated'); flag_ok = 0; end
if find(true.A_forward < repmat(bounds.a_forward_min,1,ns)), warning('Lower bounds for A violated'); flag_ok = 0; 
  display(sprintf('### Lowest A value: %4f ###',min(min(true.A_forward))));
end
if find(true.A_forward > repmat(bounds.a_forward_max,1,ns)), warning('Upper bounds for A violated'); flag_ok = 0; end


%if flag_ok == 0;     
%  warning('Artificial "true" solution violates bounds');
%end

if find(true.A_forward < 0), error('Negative driving force encountered'); end

% check whether solutions yield the right fluxes
% nn = network; 
% for it = 1:ns,
%   nn.kinetics.u = true.E(:,it);
%   [true.V(:,it), network_velocities(exp(true.X(:,it)),nn)]
% end

kinetics = network.kinetics;
[kinetics.Kcatf, kinetics.Kcatr] = modular_KV_Keq_to_kcat(network.N,network.kinetics);

% ------------------------------------------------------------
% data structure 'kinetic_data'

kinetic_data = kinetics_to_kinetic_data(network);

% ------------------------------------------------------------
% data structure 'state_data'

state_data = data_to_state_data(data);

% ------------------------------------------------------------

function x = gaussian_random_matrix(xmean,xstd,xmin,xmax,ni,nj)
  
x = xmean + xstd * randn(ni,nj);
x(x<xmin) = xmin;
x(x>xmax) = xmax;
