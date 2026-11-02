function [kinetics, prior, bounds, data, true_model, kinetic_data, state_data] = mb_generate_artificial_data(network, mb_options, q_info, c_init, conc_min, conc_max, ind_e_off, Keq_init, v_init)

% [kinetics, prior, bounds, data, true_model, kinetic_data, state_data] = mb_generate_artificial_data(network, mb_options, q_info, c_init, conc_min, conc_max, ind_e_off, Keq_init, v_init)
%
% Generate an artificial data set for a given metabolic network structure
% 
%   network: network structure
%   mb_options  struct containing options
%   mb_options.ns:      number of steady states
%   q_info       struct describing the dependencies between model variables
%   c_init       (optional) matrix; initial guess of (non-logarithmic) metabolite concentrations; 
%                           also used for generating the artificial data 
%   conc_min, conc_max vectors with minimal and maximal metabolite concentrations
%   ind_e_off    indices of enzymes to be shut off
%   Keq_init     Keq vector to be used for artificial data
%   v_init       (optional) reference flux vector to be used for artificial data
% 
%   prior        struct describing priors in the optimality problem (details see mb_make_prior)
%   bounds       struct describing bounds in the optimality problem (details see mb_make_bounds)
%   data         struct describing data used in the optimality problem
%   true_model         struct describing true model variables (optional; only for models with artificial data)

eval(default('c_init','[]', 'conc_min','[]','conc_max','[]','ind_e_off','[]','Keq_init','[]','v_init','[]'));

% --------------------------------------------------------------
%% Global variables for function modular_velocities (the Metabolic Network Toolbox currently does not read them)
global global_structure_matrices 
global_structure_matrices = 1;
global Mplus Mminus Wplus Wminus nm nr ind_M ind_Wp ind_Wm
N = network.N; W = network.regulation_matrix; ind_ext = find(network.external); h = ones(size(N,2),1);
[Mplus, Mminus, Wplus, Wminus, nm, nr, N_int,ind_M,ind_Wp,ind_Wm] = make_structure_matrices(N,W,ind_ext,h);
% --------------------------------------------------------------

[nm,nr] = size(network.N);

if ~strcmp(mb_options.ecm_score,'emc4cm'), error('wrong rate law'); end 

switch mb_options.ecm_score
  case 'emc4cm',
    network.kinetics = set_kinetics(network,'cs');
  otherwise
    error('unsupported rate laws');
end

% --------------------------------------------------------------
% kinetic constants

rng(mb_options.random_seed,'twister');

if isempty(c_init),
  c_init = ones(length(network.metabolites),mb_options.ns);
end

%% for simplicity: randomly choose std chem pot
kinetics         = network.kinetics;
mu0              = mb_options.quantities.mu0.std * randn(nm,1);
if length(Keq_init),
  display('Using predefined Keq values for artificial data');
  kinetics.Keq = Keq_init;
else
  kinetics.Keq = exp(- network.N' * mu0/RT); 
end

kinetics.KV      = exp(gaussian_random_matrix(mb_options.quantities.KV.mean_ln,mb_options.quantities.KV.std_ln,log(1.1*mb_options.quantities.KV.min),log(0.9*mb_options.quantities.KV.max),nr,1));
kinetics.KM      = double(kinetics.KM~=0) .* exp(gaussian_random_matrix(mb_options.quantities.KM.mean_ln,mb_options.quantities.KM.std_ln,log(1.1*mb_options.quantities.KM.min),log(0.9*mb_options.quantities.KM.max),nr,nm));
[kinetics.Kcatf, kinetics.Kcatr] = modular_KV_Keq_to_kcat(network.N, kinetics);

kinetics.c       = c_init;
kinetics.u       = kinetics.u .* exp(0.5*randn(nr,1));
network.kinetics = kinetics;

% check haldane relationships for sampled kinetic constants
% modular_rate_law_haldane(network)

% -------------------------------------------------------------

% [qind, qall, M_qind_to_qall, M_qall_to_qind] = mb_analyse_network(network);
% NOTE THAT THE FUNCTION NEEDS TO BE FIXED
 
% -------------------------------------------------------------
% Initialise other variables

[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);


% --------------------------------------------------------------
% Bounds

% FIX: rather get these numbers from parameter prior table

bounds = mb_make_bounds(network, q_info, mb_options, conc_min, conc_max);

% --------------------------------------------------------------
% Priors (kinetic constants)

prior = mb_make_prior(network,q_info,mb_options);

% -------------------------------------------------------------
% Generate "true" data points based on existing data set;
% use feasible state, randomise enzyme and ext. metab levels (multiplicative random variations),
% and find new steady states

clear true_model

ns = mb_options.ns;

true_model.q    = mb_kinetics_to_q(network, mb_options, q_info);
true_model.qall = mb_q_to_qall(true_model.q, q_info);

% --------------------------------------------------------------
% metabolic states

rng(mb_options.random_seed,'twister');

% compute one feasible steady state ("reference" )
e_ref = mb_options.metabolic_prior_e_geom_mean * [mb_options.metabolic_artificial_e_geom_std .^ randn(nr,1)];
e_ref(ind_e_off,:) = 10^-15;

network.kinetics.u = e_ref;
kinetics.u         = e_ref;

if length(v_init),
  display('mb_generate_artificial_data: adjusting reference state to given fluxes and concentrations')
  % VERSION 1: USE original fluxes and metabolite concentrations directly and determine necessary enzyme levels
  % -> later there are still flux reversals
  myv = network_velocities(c_init,network,kinetics);
  c_ref = c_init;
  e_ref = e_ref .* v_init ./ myv;
  v_ref = v_init;
  %check forces and fluxes: 
  %[log(Keq_init)-network.N'*log(c_init), v_init, myv, v_ref]

  % % VERSION 2: USE original fluxes and run ECM to obtain metabolite + enzyme concentrations for them
  % % -> with suitable (small) variation between states there are no flux reversals, but met and enz conc exceeed upper limits
  % % -> not feasible in MB! 
  % ecm_options            = ecm_default_options(network, 'E. coli central carbon metabolism');
  % ecm_options.ecm_scores = {'emc4cm'};
  % ecm_options            = ecm_update_options(network, ecm_options);
  % [c_struct, u_struct] = ecm_enzyme_cost_minimization(network, network.kinetics, v_init, ecm_options);
  % c_ref = c_struct.emc4cm;
  % e_ref = u_struct.emc4cm;

else
  [c_ref, v_ref] = network_steady_state(network, network.kinetics.c, 100000);
end

x_ref  = log(c_ref);

% check correctness and driving forces
% pp   = mb_make_pp(network);
% pp.v = v_ref;
% [~, e_ref] = ecm_get_score(mb_options.ecm_score,x_ref,pp);
% [network.kinetics.u, e_ref]
% theta = driving_force(c_ref,v_ref,network)
% pm(RT * theta, network.actions)

% perturb the reference state to obtain several steady states ("samples")

X_init =     repmat(x_ref,1,ns)      + log(mb_options.metabolic_artificial_c_geom_std) * randn(nm,ns);
true_model.E = exp(repmat(log(e_ref),1,ns) + log(mb_options.metabolic_artificial_e_geom_std) * randn(nr,ns));
true_model.E(e_ref==0) = 0;

for j=1:ns,
  my_network = network;
  my_network.kinetics.c = exp(X_init(:,j));
  my_network.kinetics.u = true_model.E(:,j);
  [c_ss,v_ss] = network_steady_state(my_network, my_network.kinetics.c);
  %check stationarity: network.N(network.external==0,:)*v_ss
  true_model.X(:,j) = log(c_ss);
  true_model.V(:,j) = v_ss;
  true_model.A(:,j) = RT * [log(network.kinetics.Keq) - network.N' * true_model.X(:,j)]; % kJ/mol
  true_model.A_forward(:,j) = sign(sign(v_ss)+0.5) .* true_model.A(:,j);
end

if find(true_model.A .* true_model.V <0), error('infeasible flux direction'); end

%  log(network.kinetics.Keq)-network.N'*log(c_ref)

true_model.kinetics = network.kinetics;


% --------------------------------------------------------------
% Data matrices

for it = 1:ns,
  data.samples{it,1} = ['state_' num2str(it)];
end

data.V.mean = true_model.V;
data.X.mean = true_model.X;
data.lnE.mean = log(true_model.E);
%data.lnE.mean(true_model.E==0) = 10^-15;

data.X.std   = log(mb_options.data_C_geom_std)    * ones(size(data.X.mean));
data.lnE.std = log(mb_options.data_E_geom_std)    * ones(size(data.lnE.mean));
data.V.std   = [mb_options.data_V_geom_std-1]  * abs(data.V.mean) + 0.01 * max(max(abs(data.V.mean)));

% noise in artificial data (NOTE: the resulting "states" may be thermo-physiologically infeasible!)

if mb_options.use_artificial_noise,
  if mb_options.verbose,
    display('Using artificial metabolic data with noise');
  end
  ind_zero = find(data.V.mean==0);
  data.V.mean = data.V.mean + data.V.std .* randn(nr,ns);
  data.V.mean(ind_zero)=0;
  ind_wrong_sign = find(data.V.mean .* true_model.V <0);
  data.V.mean(ind_wrong_sign) = -data.V.mean(ind_wrong_sign);
  data.X.mean = data.X.mean + data.X.std .* randn(nm,ns);
  data.lnE.mean = data.lnE.mean + data.lnE.std .* randn(nr,ns);
else
  if mb_options.verbose,
    display('Using artificial metabolic data without noise');
  end
end

data.qall.mean = mb_q_to_qall(true_model.q,q_info);
data.qall.std  = log(mb_options.data_kin_geom_std) * ones(size(data.qall.mean));

if mb_options.use_kinetic_data_noise,
  if mb_options.verbose,
    display('Generating artificial kinetic data with noise');
  end
  data.qall.mean = data.qall.mean + data.qall.std .* randn(size(data.qall.mean));
else
  if mb_options.verbose,
    display('Generating artificial kinetic data without noise');
  end
end

data.kinetics_median = mb_qall_to_kinetics(data.qall.mean,network,mb_options,q_info);

% --------------------------------------------------------------
% check whether true values satisfy bounds

% flag_ok = 1;
if find(true_model.q < bounds.q_min),       error('Artificial data violate Lower bounds for q - MB will not work'); end
if find(true_model.q > bounds.q_max),               error('Artificial data violate Upper bounds for q - MB will not work'); end
if find(true_model.X < repmat(bounds.x_min,1,ns)),  error('Artificial data violate Lower bounds for X - MB will not work'); end
if find(true_model.X > repmat(bounds.x_max,1,ns)),  error('Artificial data violate Upper bounds for X - MB will not work'); end
if find(true_model.E < repmat(bounds.e_min,1,ns)),  error('Artificial data violate Lower bounds for E - MB will not work'); end
if find(true_model.E > repmat(bounds.e_max,1,ns)), 
  factor = max(max(true_model.E./repmat(bounds.e_max,1,ns))) 
  error('Artificial data violate Upper bounds for E - MB will not work'); end
if find(true_model.A_forward < repmat(bounds.a_forward_min,1,ns)), error('Sampled states violate lower bound for reaction affinity A - please consider using a different initial concentration profile'); flag_ok = 0; 
  % display(sprintf('### Lowest A value: %4f ###',min(min(true_model.A_forward))));
end
if find(true_model.A_forward > repmat(bounds.a_forward_max,1,ns)), warning('Sampled states violate upper bound for reaction affinity A - please consider using a different initial concentration profile'); flag_ok = 0; end


%if flag_ok == 0;     
%  warning('Artificial "true" solution violates bounds');
%end

if find(true_model.A_forward < 0), error('Negative driving force encountered'); end

% check whether solutions yield the right fluxes
% nn = network; 
% for it = 1:ns,
%   nn.kinetics.u = true_model.E(:,it);
%   [true_model.V(:,it), network_velocities(exp(true.X(:,it)),nn)]
% end

% ------------------------------------------------------------
% data structure 'kinetic_data'

kinetic_data = kinetics_to_kinetic_data(network);

% ------------------------------------------------------------
% data structure 'state_data'

state_data = data_to_state_data(data);

%% Clear global variables
clearvars -global global_structure_matrices Mplus Mminus Wplus Wminus nm nr ind_M ind_Wp ind_Wm


% ------------------------------------------------------------

function x = gaussian_random_matrix(xmean,xstd,xmin,xmax,ni,nj)
  
x = xmean + xstd * randn(ni,nj);
x(x<xmin) = xmin;
x(x>xmax) = xmax;
