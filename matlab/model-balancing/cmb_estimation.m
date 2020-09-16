function [optimal, gradient_V, init, time] = cmb_estimation(network, q_info, bounds, data, prior, preposterior, pp, V, cmb_options)

% [optimal, gradient_V, init, time] = cmb_estimation(network, q_info, bounds, data, prior, preposterior, pp, V, cmb_options)
%
% gradient_V Gradient of enzyme posterior loss with respect to fluxes (currently not computed)
%  
% --------------------------------------------------------
% Initial guesses for X and q (where preposterior q means: using prior (and possibly data) for q only 
% 
% X_init = argmin_X sum(sum( [ [X - preposterior.X.mean] ./ preposterior.X.std ] .^2));  WITHIN allowed polytopes!
% q_init = argmin_X sum(sum( [ [q - preposterior.q.mean] ./ preposterior.q.std     ] .^2);   WITHIN allowed ranges!
% --------------------------------------------------------

% --------------------------------------------------------------
%% Use global variables to speed up function modular_velocities
global global_structure_matrices
global Mplus Mminus Wplus Wminus nm nr ind_M ind_Wp ind_Wm
global LP_info % variable used in cmb_log_posterior
% --------------------------------------------------------------

N = network.N; W = network.regulation_matrix; ind_ext = find(network.external); h = network.kinetics.h;
[Mplus, Mminus, Wplus, Wminus, nm, nr, N_int,ind_M,ind_Wp,ind_Wm] = make_structure_matrices(N,W,ind_ext,h);

% --------------------------------------------------------------

% --------------------------------------------------------
% Initialise some variables

tstart = tic;

epsilon = 10^-15;

[nr,nm,nx,KM_indices,KA_indices,KI_indices,nKM,nKA,nKI] = network_numbers(network);
[~,ns] = size(data.X.mean);
nq     = length(prior.q.mean);

Aforward_min = cmb_options.quantities.Aforward.min * ones(nr,ns);


% --------------------------------------------------
% Bounds and preposterior distributions for y vector

y_bound_min            = cmb_qX_to_y(bounds.q_min, repmat(bounds.x_min,1,ns), nm,ns);
y_bound_max            = cmb_qX_to_y(bounds.q_max, repmat(bounds.x_max,1,ns), nm,ns);
y_preposterior_mean    = cmb_qX_to_y(preposterior.q.mean, preposterior.X.mean, nm, ns);
y_preposterior_cov_inv = preposterior.q.cov_inv; 
for it = 1:ns,
  y_preposterior_cov_inv = matrix_add_block(y_preposterior_cov_inv, diag(1./preposterior.X.std(:,it).^2));
end

% --------------------------------------------------
% Build matrix describing the flux sign constraints for y. 
%
% consists of vertically stacked blocks; each block describes the flux sign constraints for one state
% in each block: put subblocks blocks together horizontally: first X then q, as in cmb_qX_to_y

block_X = -network.N';
block_q = q_info.M_q_to_qall(q_info.qall.index.Keq,:);
y_ineq_A = [];
for it = 1:ns,
  y_ineq_A = matrix_add_block(y_ineq_A, block_X);
end

y_ineq_A = [y_ineq_A, repmat(block_q, ns, 1)];
y_ineq_A = -diag(sign(V(:))) * y_ineq_A;
y_ineq_b = -[cmb_options.quantities.Aforward.min * ones(nr * ns,1)];

ind_finite_flux = find(V(:)~=0);
y_ineq_A = y_ineq_A(ind_finite_flux,:);
y_ineq_b = y_ineq_b(ind_finite_flux,:);


% --------------------------------------------------
% Set initial values

switch cmb_options.initial_values_variant,
  case {'given_values','true_values','random'},

    y_init = cmb_qX_to_y(cmb_options.init.q, cmb_options.init.X, nm,ns);
    init   = cmb_options.init;
  
  otherwise,
    %% if y_preposterior_mean satisfies all constraints, use it
    if 0 == sum(y_preposterior_mean < y_bound_min) ...
          + sum(y_preposterior_mean > y_bound_max)  ...
          + sum(y_ineq_A * y_preposterior_mean > y_ineq_b),
      y_init = y_preposterior_mean;
    else
      %% find preposterior mode under constraints
      %% FIX ME, use cplexqp (to account for constraints)
      opt_quadprog = struct('Algorithm','interior-point-convex','MaxFunEvals',10^10,'MaxIter',10^10,'TolX',10^-5,'Display','off','OptimalityTolerance', 10^-5,'StepTolerance', 10^-10);
      [y_init,~,err] = quadprog(y_preposterior_cov_inv, -y_preposterior_cov_inv * y_preposterior_mean, y_ineq_A, y_ineq_b-epsilon,[],[], y_bound_min, y_bound_max,[],opt_quadprog);
      if err<0, error(sprintf('Error in computing initial values: quadprog error flag %d',err)); end
    end
    [init.q, init.X] = cmb_y_to_qX(y_init, nm, ns);

end

% ------------------------------------------------------
% Check whether initial values satify all constraints

init_feasible = 1;

if sum(y_init < y_bound_min) + sum(y_init > y_bound_max) ~= 0,
  warning('Initial point violates some bounds');
  init_feasible = 0;
end

if length(y_ineq_b),
  if sum(y_ineq_A * y_init > y_ineq_b-epsilon) ~=0,
    warning('Initial point violates some constraints; please change the initial state (or possibly Aforward_min)');
    occuring_values_and_constraints = [y_ineq_A * y_init, y_ineq_b]
    init_feasible = 0;
  end
end

if init_feasible,
  if cmb_options.verbose,
    display('Initial point respects all bounds and constraints');
  end
end

% ------------------------------------------------------
% Set tighter bounds, derived from initial objective value and the conditions

% 0.5 * [new_X_bound -preposterior.X.mean].^2 ./ preposterior.X.std.^2 < f_init
% 0.5 * [new_q_bound -preposterior.q.mean].^2 ./ preposterior.q.std.^2 < f_init

LP_info.y_ineq_A = y_ineq_A;
LP_info.y_ineq_b = y_ineq_b;
LP_info.epsilon  = epsilon;

if init_feasible,
  f_init = cmb_log_posterior(y_init,pp,preposterior,V,cmb_options,q_info);
  preposterior.q.std = sqrt(diag(inv(preposterior.q.cov_inv)));
  new_X_min       = preposterior.X.mean - sqrt(2 * abs(f_init)) * preposterior.X.std;
  new_X_max       = preposterior.X.mean + sqrt(2 * abs(f_init)) * preposterior.X.std;
  new_q_min       = preposterior.q.mean - sqrt(2 * abs(f_init)) * preposterior.q.std;
  new_q_max       = preposterior.q.mean + sqrt(2 * abs(f_init)) * preposterior.q.std;
  y_bound_min_new = cmb_qX_to_y(bounds.q_min, repmat(bounds.x_min,1,ns), nm,ns);
  y_bound_max_new = cmb_qX_to_y(bounds.q_max, repmat(bounds.x_max,1,ns), nm,ns); 
  y_bound_min     = max(y_bound_min,y_bound_min_new);
  y_bound_max     = min(y_bound_max,y_bound_max_new);
end

% --------------------------------------------------------
% Optimisation

% TEST: show log posterior and log posterior gradient at the q,X-preposterior mean
% y_mean = cmb_qX_to_y(preposterior.q.mean,preposterior.X.mean,nm,ns);
% cmb_log_posterior(y_mean,pp,preposterior,V,cmb_options,q_info);

switch cmb_options.score,
  case 'log_neg_log_posterior',
  otherwise
    error('optimisation of neg_log_posterior not implemented')
end

if cmb_options.use_gradient,
  opt = optimoptions('fmincon','MaxFunEvals',10^15,'MaxIter',10^15,'TolX',10^-5,'Display',...
                     cmb_options.optim_display,'Algorithm','interior-point','SpecifyObjectiveGradient',true);
  [y_opt,~,err] = fmincon(@(y) cmb_objective(y,pp,preposterior,V,cmb_options,q_info),y_init,y_ineq_A,y_ineq_b-epsilon,[],[],y_bound_min,y_bound_max,[],opt);
else
  opt = optimoptions('fmincon','MaxFunEvals',10^15,'MaxIter',10^15,'TolX',10^-5,'Display',...
                     cmb_options.optim_display,'Algorithm','interior-point','SpecifyObjectiveGradient',false);
  [y_opt,~,err] = fmincon(@(y) log(-cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,cmb_options.verbose)),y_init,y_ineq_A,y_ineq_b-epsilon,[],[],y_bound_min,y_bound_max,[],opt);
end

if err<0,
  error(sprintf('No feasible solution found. Error flag %d',err));
end

%CHECK: show function value of initial guess and end result
f_init = cmb_log_posterior(y_init,pp,preposterior,V,cmb_options,q_info);
f_opt  = cmb_log_posterior(y_opt,pp,preposterior,V,cmb_options,q_info);
if f_init>f_opt,
  if init_feasible,
    display('Initial guess was better than optimised value ("THANK YOU" matlab); using initial value');
    y_opt = y_init;
  end
end

[q_opt,X_opt] = cmb_y_to_qX(y_opt,nm,ns);

optimal = cmb_qX_to_variables(q_opt,X_opt,V,network,cmb_options,q_info,cmb_options.ecm_score,pp);

% --------------------------------------------------------
% Gradient of enzyme posterior loss with respect to fluxes

% CHECK ME AGAIN
gradient_V = [[optimal.E - preposterior.E.mean] ./ preposterior.E.std] .* [optimal.E ./ V];
gradient_V = nan; 

time = toc(tstart);

display(sprintf('Calculation time: %3.2f s',time));

