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
% 
% choice of initial values depends on cmb_options.initial_values_variant:
%  'average_sample', 
%  'preposterior_mode', 
%  'random'
%  'true_values'
%  'given_values'

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

tstart_cmb_estimation = tic;

epsilon = 10^-9;

[nr,nm,nx,KM_indices,KA_indices,KI_indices,nKM,nKA,nKI] = network_numbers(network);
[~,ns] = size(data.X.mean);
nq     = length(prior.q.mean);

Aforward_min = cmb_options.quantities.Aforward.min * ones(nr,ns);


% --------------------------------------------------
% Bounds and preposterior distributions for y vector

y_bound_min = cmb_qX_to_y(bounds.q_min, repmat(bounds.x_min,1,ns), nm,ns);
y_bound_max = cmb_qX_to_y(bounds.q_max, repmat(bounds.x_max,1,ns), nm,ns);

y_preposterior_mean    = cmb_qX_to_y(preposterior.q.mean, preposterior.X.mean, nm, ns);
y_preposterior_cov_inv = preposterior.q.cov_inv; 
for it = 1:ns,
  y_preposterior_cov_inv = matrix_add_block(y_preposterior_cov_inv, diag(1./preposterior.X.std(:,it).^2));
end

% --------------------------------------------------
% Build matrix describing the flux sign constraints for y. 
%
% consists of vertically stacked blocks; each block describes the flux sign constraints for one state
% in each block: put sub-blocks together horizontally: first X then q, as in cmb_qX_to_y

block_X = -network.N';
block_q = q_info.M_q_to_qall(q_info.qall.index.Keq,:);
y_ineq_A = [];

for it = 1:ns,
  y_ineq_A = matrix_add_block(y_ineq_A, block_X);
end

y_ineq_A = [y_ineq_A, repmat(block_q, ns, 1)];
y_ineq_A = -diag(sign(V(:))) * y_ineq_A;  
% note that this correctly implements the fact that in reactions with v==0 there is no constraint
y_ineq_b = -[cmb_options.quantities.Aforward.min * ones(nr * ns,1)];

ind_finite_flux = find(V(:)~=0);
y_ineq_A = y_ineq_A(ind_finite_flux,:);
y_ineq_b = y_ineq_b(ind_finite_flux,:);

% constraints for bounds on dependent kinetic constants

y_ineq_bounds_A = y_ineq_A;
y_ineq_bounds_b = y_ineq_b;

% cmb_options.use_bounds == 0, the constraints from bounds are simply ignored

nq       = q_info.q.number;
M_y_to_q = [zeros(nq,nm*ns), eye(nq)];
YY_A     = [-q_info.M_q_to_qall * M_y_to_q; ...
             q_info.M_q_to_qall * M_y_to_q];
YY_b     = [-bounds.q_all_min; ...
             bounds.q_all_max];

y_ineq_bounds_A = [y_ineq_bounds_A; YY_A];
y_ineq_bounds_b = [y_ineq_bounds_b; YY_b];

if cmb_options.use_bounds,
  y_ineq_A = y_ineq_bounds_A;
  y_ineq_b = y_ineq_bounds_b;
end

% --------------------------------------------------
% Set initial values

switch cmb_options.initial_values_variant,
  
  case 'polytope center',
    %% find initial point using linprog
    nn = 5;
    y_init = find_polytope_centre_randn([],[],y_ineq_bounds_A,y_ineq_bounds_b,y_bound_min, y_bound_max, epsilon, nn);

    %% %%   % check validity of the solution
    %% %y_bound_min-y_init
    %% find(y_init<y_bound_min)
    %% %y_init-y_bound_max
    %% find(y_init>y_bound_max)
    %% %max(y_ineq_bounds_A * y_init - y_ineq_bounds_b)
    %% find(y_ineq_bounds_A * y_init > y_ineq_bounds_b)

    %% clear y_init_list
    %% for it =1:nn;
    %%   
    %%   f = randn(size(y_bound_min));
    %%   opt = optimoptions('linprog','Display','none','Algorithm','interior-point','ConstraintTolerance',10^-10);
    %%   [y_init_list(:,it), ~,exitflag]   = linprog(f,y_ineq_bounds_A,y_ineq_bounds_b-epsilon,[],[],y_bound_min+epsilon,y_bound_max-epsilon,[],opt);
    %%   [y_init_list(:,it+nn),~,exitflag] = linprog(-f,y_ineq_bounds_A,y_ineq_bounds_b-epsilon,[],[],y_bound_min+epsilon,y_bound_max-epsilon,[],opt);
    %% 
    %%   % check validity of the solution
    %%   % y_bound_min-y_init_list(:,it)
    %%   % y_init_list(:,it)-y_bound_max
    %%   % y_bound_min-y_init_list(:,it+nn)
    %%   % y_init_list(:,it+nn)-y_bound_max
    %%   % max(y_ineq_bounds_A * y_init_list(:,it)- y_ineq_bounds_b)
    %%   % max(y_ineq_bounds_A * y_init_list(:,it+nn) - y_ineq_bounds_b)
    %%   
    %%   my_init_q = cmb_y_to_qX(y_init_list(:,it+nn), nm, ns);
    %%   if sum(q_info.M_q_to_qall * my_init_q < bounds.q_all_min),
    %%     q_info.M_q_to_qall * my_init_q - bounds.q_all_min
    %%     error('Infeasible initial point'); end
    %%   if sum(q_info.M_q_to_qall * my_init_q > bounds.q_all_max), 
    %%     q_info.M_q_to_qall * my_init_q - bounds.q_all_max
    %%     error('Infeasible initial point'); end
    %% 
    %% end
    %% y_init = mean(y_init_list,2);

    if prod(y_ineq_bounds_A * y_init < y_ineq_bounds_b) * prod(y_init >= y_bound_min) * prod(y_init <= y_bound_max) ==0,
      error('Infeasible solution');
    end
    [init.q, init.X] = cmb_y_to_qX(y_init, nm, ns);
  
  case 'flat_objective',
    %% find initial point using fmincon with a flat objective function
    opt = optimoptions('fmincon','MaxFunEvals',10^15,'MaxIter',10^15,'TolX',10^-5,'Display',...
                       cmb_options.optim_display,'Algorithm','interior-point','SpecifyObjectiveGradient',false,'ConstraintTolerance',10^-10);
    [y_init,~,err] = fmincon(@(y) 1, y_init, y_ineq_bounds_A, y_bounds_ineq_b-epsilon, [],[],y_bound_min,y_bound_max,[],opt);
    [init.q, init.X] = cmb_y_to_qX(y_init, nm, ns);
    
  case {'given_values','true_values','random'},
    init   = cmb_options.init;
    y_init = cmb_qX_to_y(init.q, init.X, nm,ns);

  otherwise,
    % if y_preposterior_mean satisfies all constraints, use it
    if 0 == sum(y_preposterior_mean < y_bound_min) ...
          + sum(y_preposterior_mean > y_bound_max)  ...
          + sum(y_ineq_bounds_A * y_preposterior_mean > y_ineq_bounds_b),
      y_init = y_preposterior_mean;
    else
      %% find preposterior mode under constraints
      %% in the future, use cplexqp (to account for constraints)
      opt_quadprog = struct('Algorithm','interior-point-convex','MaxFunEvals',10^10,'MaxIter',10^10,'TolX',10^-5,'Display','off','OptimalityTolerance', 10^-5,'StepTolerance', 10^-10);
      [y_init,~,err] = quadprog(y_preposterior_cov_inv, -y_preposterior_cov_inv * y_preposterior_mean, y_bounds_ineq_A, y_bounds_ineq_b-epsilon,[],[], y_bound_min, y_bound_max,[],opt_quadprog);
      if err<0, error(sprintf('Error in computing initial values: quadprog error flag %d',err)); end
    end
    [init.q, init.X] = cmb_y_to_qX(y_init, nm, ns);

end

% -----------------------------------------------
% Check feasibility of initial point (check explicitly with init.q and init.X)

% Kinetic constants, bounds 

if cmb_options.use_bounds,
  if sum(q_info.M_q_to_qall * init.q < bounds.q_all_min), 
    q_info.M_q_to_qall * init.q - bounds.q_all_min
    error('Infeasible initial point'); 
  end
  if sum(q_info.M_q_to_qall * init.q > bounds.q_all_max), 
    q_info.M_q_to_qall * init.q - bounds.q_all_max
    error('Infeasible initial point'); 
  end
end

% Metabolite data, bounds 

if cmb_options.use_bounds,
  if sum(sum(init.X < repmat(bounds.x_min,1,ns))), error('Infeasible initial point'); end
  if sum(sum(init.X > repmat(bounds.x_max,1,ns))), error('Infeasible initial point'); end
end

% Thermodynamic forces, directions

my_ln_Keq       = init.q(q_info.q.index.Keq);
my_theta_matrix = repmat(my_ln_Keq,1,ns) - network.N' * init.X;

if sum(sum(my_theta_matrix .* data.V.mean < 0)),          error('Infeasible initial point'); end
if sum(sum([my_theta_matrix == 0] .* [data.V.mean ~=0])), error('Infeasible initial point'); end


% ------------------------------------------------------
% Check whether initial values satify all constraints (check with vector y_init)

init_feasible = 1;

if cmb_options.use_bounds,
  if sum(y_init < y_bound_min) + sum(y_init > y_bound_max) ~= 0,
    error('Initial point violates some bounds');
    init_feasible = 0;
  end
end

if length(y_ineq_b),
  if sum(y_ineq_A * y_init > y_ineq_b) ~=0, % -epsilon
    error('Initial point violates some constraints; please change the initial state (or possibly Aforward_min)');
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

if cmb_options.use_bounds,
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
end


% --------------------------------------------------------
% Optimisation

% TEST: show log posterior and log posterior gradient at the q,X-preposterior mean
% y_mean = cmb_qX_to_y(preposterior.q.mean,preposterior.X.mean,nm,ns);
% cmb_log_posterior(y_mean,pp,preposterior,V,cmb_options,q_info);

if cmb_options.use_bounds == 0, 
  y_bound_min = -inf * ones(size(y_bound_min));
  y_bound_max = inf * ones(size(y_bound_min));
end

opt = optimoptions('fmincon','MaxFunEvals',10^15,'MaxIter',10^15,'TolX',10^-5,'Display',...
                   cmb_options.optim_display,'Algorithm','interior-point','SpecifyObjectiveGradient',...
                   false,'ConstraintTolerance',10^-10);

if cmb_options.use_gradient,
  opt.SpecifyObjectiveGradient = true;
end

[y_opt,~,err] = fmincon(@(y) cmb_objective(y,pp,preposterior,V,cmb_options,q_info,cmb_options.verbose),y_init,y_ineq_A, y_ineq_b-epsilon, [],[],y_bound_min,y_bound_max,[],opt);

% % check validity of the solution (violation of inequality constraints
% err
% epsilon
% %y_bound_min-y_opt % must be negative!
% %y_opt-y_bound_max % must be negative!
% y_ineq_A * y_opt - y_ineq_b % must be negative!

if err<0,
  error(sprintf('No feasible solution found. Error flag %d',err));
end

%CHECK: show function value of initial guess and end result
f_init = cmb_log_posterior(y_init,pp,preposterior,V,cmb_options,q_info);
f_opt  = cmb_log_posterior(y_opt,pp,preposterior,V,cmb_options,q_info);
if f_init>f_opt,
  if init_feasible,
    display('Optimised value is worse than initial guess; keeping initial value');
    y_opt = y_init;
  end
end

[q_opt,X_opt] = cmb_y_to_qX(y_opt,nm,ns);

optimal = cmb_qX_to_variables(q_opt,X_opt,V,network,cmb_options,q_info,cmb_options.ecm_score,pp);

if cmb_options.use_bounds,
  if sum(q_info.M_q_to_qall * q_opt < bounds.q_all_min), 
    q_info.M_q_to_qall * q_opt - bounds.q_all_min
    error('Infeasible solution point'); 
  end
  if sum(q_info.M_q_to_qall * q_opt > bounds.q_all_max), 
    q_info.M_q_to_qall * q_opt - bounds.q_all_max
    error('Infeasible solution point'); 
  end
end


% --------------------------------------------------------
% Gradient of enzyme posterior loss with respect to fluxes

% CHECK ME AGAIN
%gradient_V = [[optimal.E - preposterior.E.mean] ./ preposterior.E.std] .* [optimal.E ./ V];
gradient_V = [[log(optimal.E) - preposterior.lnE.mean] ./ preposterior.lnE.std] .* [optimal.E ./ V];
gradient_V = nan; 

time = toc(tstart_cmb_estimation);
display(sprintf('\nCalculation time: %3.2f s',time));

% --------------------------------------------------------------
%% Clear global variables
clearvars -global global_structure_matrices Mplus Mminus Wplus Wminus nm nr ind_M ind_Wp ind_Wm 
%LP_info
% --------------------------------------------------------------
