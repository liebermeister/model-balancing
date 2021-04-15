function [optimal, calculation_time, gradient_V, init, cmb_options, V, kapp_max, preposterior, pp] = model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, true, init)

% [optimal, calculation_time, gradient_V, init, cmb_options, V, kapp_max, preposterior, pp] = model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, true, init)
%
% Convex model balancing
% For usage examples, see the demo files 'demo_cmb_artificial_data.m' and 'demo_cmb_experimental_data.m'
%
% Input variables:
%   filenames    struct, (model- and run-specific) filenames for input and output files
%     .graphics_dir optional; output directory for saving graphics
%     .report_txt   optional; output report file (.txt)
%     .results_mat  optional;  output results file (.mat)
%   cmb_options  struct, options 
%   network      struct, metabolic network (format as in Metabolic Network Toolbox)
%   q_info       struct, the dependencies between model variables
%   prior        struct, priors in the optimality problem (details see cmb_make_prior)
%   bounds       struct, bounds in the optimality problem (details see cmb_make_bounds)
%   data         struct, data used in the optimality problem
%                  flux maxtrix: data.V.mean
%   true         struct, true model variables (optional; only for models with artificial data)
%   init         struct, initial values
%
% Output variables:
%   optimal          struct containing the model variables obtained from optimisation
%   calculation_time calculation time in seconds
%   gradient_V       gradient of enzyme posterior loss w.r.t. fluxes (currently not computed)
%   init             struct containing the initial model variables (fields init.q and init.X)
%   cmb_options      struct containing options (updated by the function)
%   V                flux data 
%   kapp_max         struct containing kapp_max values
%   preposterior     struct containing the preposterior distribution
%   pp               struct containing some data used during calculation (see cmb_make_pp(network))
%
% For generating the input data, see 'cmb_model_and_data' and 'cmb_model_artificial_data'

if cmb_options.display,
  display(sprintf('\n---------------'));
  display(sprintf('Model balancing'));
  display(sprintf('---------------\n'));
end
  
filenames_default = struct('graphics_dir',[],'report_txt',[],'results_mat',[]);
filenames = join_struct(filenames_default,filenames);
  
tstart_model_balancing = tic;

  % --------------------------------------------------------------
%% Set global variables to speed up function modular_velocities
global global_structure_matrices
global_structure_matrices = 1;
global Mplus Mminus Wplus Wminus nm nr ind_M ind_Wp ind_Wm
global LP_info % variable used in cmb_log_posterior

% --------------------------------------------------------------

eval(default('true','[]'));

ns = size(data.X.mean,2);
cmb_options.ns = ns;

% -----------------------------------------------

switch cmb_options.enzyme_score_type,
  case 'quadratic'
    % normal quadratic likelihood term: with this option, MB is not guaranteed to be convex!
    %display('model_balancing.m:');
    display('Using quadratic formula for enzyme posterior score.');
    %display('  The optimality problem may be non-convex');
  case 'monotonic'
    %display('model_balancing.m:')
    display('Using convex formula for enzyme posterior score - enzyme levels may be underestimated.')
    %display('  The optimality problem is convex, but enzyme levels may be underestimated!');
  case 'interpolated'
    %display(sprintf('model_balancing.m:'))
    display(sprintf('Using enzyme fit stringency alpha=%f ',cmb_options.enzyme_score_alpha))
    display('The optimality problem may be non-convex, and enzyme levels may be underestimated.');
end


% -----------------------------------------------
% Initialise some variables

[nr,nm] = network_numbers(network);

ns = size(data.X.mean,2);
nq = length(prior.q.mean);

% Prepare the preposterior distributions (individual terms)
preposterior = cmb_prepare_posterior(prior, data, cmb_options, q_info);

if sum(isnan(data.V.mean(:))),
  error('Some fluxes are unknown');
end

V = data.V.mean;

pp = cmb_make_pp(network);

if isempty(true),
  if strcmp(cmb_options.initial_values_variant,'true_values'),
    warning('No true values given; using method "preposterior_mode" for choosing initial values');
    cmb_options.initial_values_variant = 'preposterior_mode';
  end
end

conc_min = bounds.conc_min;
conc_max = bounds.conc_max;

% check whether fluxes are thermodynamically feasible
thermo_pb_options.c_min = conc_min;
thermo_pb_options.c_max = conc_max;
for it = 1:size(V,2),
  [~,~,~, feasible] = thermo_pb(network.N, V(:,it), thermo_pb_options, 0);
  if ~feasible,
    error('Flux distribution is thermodynamically infeasible');
  end
end
display('Flux distributions are thermo-physiologically feasible');


% -----------------------------------------------
% Initial values

switch cmb_options.initial_values_variant,
  case 'polytope center';
    display(sprintf('Using center of mass of random LP solutions as the initial point'));
    cmb_options.init = [];
  case 'flat_objective';
    display(sprintf('Using fmincon solution with flat objective (under constraints) as the initial point'));
    cmb_options.init = [];
  case 'preposterior_mode';
    display(sprintf('Using preposterior mode (under constraints) as the initial point'));
    cmb_options.init = [];
  case 'random',
    display('Using random values close to zero as the initial point; these values may violate the constraints');
    cmb_options.init.q = 0.01 * randn(nq,1);
    cmb_options.init.X = 0.01 * randn(size(data.X.mean));
  case 'true_values',
    display('Using true values as the initial point');
    cmb_options.init.q = true.q;
    cmb_options.init.X = true.X;
  case 'given_values';
    display('Using given initial values');
    if ~isfield(cmb_options,'init'), 
      error('missing options field "cmb_options.init", required by "cmb_options.initial_values_variant=given_values"'); 
    end
  case 'average_sample';
    display(sprintf('Using model-balancing random data as the initial point'));

    my_cmb_options                        = cmb_options;
    my_cmb_options.ns                     = 1;
    my_cmb_options.initial_values_variant = 'preposterior_mode'; 
    my_cmb_options.show_graphics          = 0;
    my_cmb_options.verbose                = 0;
    my_cmb_options.display                = 0;
    my_cmb_options.save_results           = 0;
    
    my_q_info        = cmb_define_parameterisation(network, my_cmb_options); 
    [~, my_prior, my_bounds, my_data] = cmb_generate_artificial_data(network, my_cmb_options, my_q_info,ones(size(network.metabolites)), conc_min, conc_max);
    
    my_data.V.mean = nanmean(data.V.mean,2);
    my_data.X.mean = nanmean(data.X.mean,2);
    my_data.E.mean = exp(nanmean(data.lnE.mean,2));
    %% taking the average over standard deviations is not ideal .. maybe change this?
    my_data.V.std  = nanmean(data.V.std,2); 
    my_data.X.std  = nanmean(data.X.std,2);
    my_data.E.std  = exp(nanmean(data.lnE.std,2));
    
    if length(true),
      my_true   = true;
      my_true.V = nanmean(true.V,2);
      my_true.X = nanmean(true.X,2);
      my_true.E = exp(nanmean(log(true.E),2));
    else
      my_true = [];
    end
    
    my_optimal = model_balancing(filenames, my_cmb_options, network, my_q_info, my_prior, my_bounds, my_data, my_true);
    cmb_options.init.q = my_optimal.q;
    cmb_options.init.X = repmat(my_optimal.X,1,cmb_options.ns);
    %display(sprintf('Initial point found\n-------------------'));
  
  otherwise 
    error('unknown option');
end

if length(cmb_options.init),
  cmb_options.init.q(isnan(cmb_options.init.q)) = 0;
  cmb_options.init.X(isnan(cmb_options.init.X)) = 0;
end


% -----------------------------------------------
% Compute the optimum state and parameter values

if cmb_options.verbose,
  if cmb_options.use_gradient,
    display('Using gradients in optimisation');
  else
    display('Not using gradients in optimisation');
  end
end


display(' '); 
display('Running optimisation ..');

[optimal, gradient_V, init] = cmb_estimation(network, q_info, bounds, data, prior, preposterior, pp, V, cmb_options);

optimal.A_forward = optimal.A .* sign(sign(V)+0.5);


% -----------------------------------------------
% Compare posterior loss scores for true values, 
% initial guess, and predicted solution
% -----------------------------------------------

if cmb_options.display,
  cmb_display_scores(network, q_info, cmb_options, pp, preposterior, init, optimal, true, V, cmb_options.verbose);
end


% -----------------------------------------------
% Further improve the solution (or try it as a check) 
% by additional optimisation rounds
% -----------------------------------------------

if cmb_options.use_safe_optimisation,
  cont = 1;
  my_cmb_options = cmb_options;
  while cont,
    my_cmb_options.init = optimal;  
    my_cmb_options.initial_values_variant = 'given_values';
    [optimal, gradient_V, my_init] = cmb_estimation(network, q_info, bounds, data, prior, preposterior, pp, V, my_cmb_options);
    cont = cmb_display_scores(network, q_info, my_cmb_options, pp, preposterior, my_init, optimal, true, V, cmb_options.verbose);
  end
end


% ---------------------------------------------------------------------
% Estimation of kcat values by maximal kapp values (as in Davidi et al)
% ---------------------------------------------------------------------

%kapp_max.forward =  max([data.V.mean ./ data.E.mean],[],2); 
%kapp_max.reverse = -min([data.V.mean ./ data.E.mean],[],2);
kapp_max.forward =  max([data.V.mean ./ exp(data.lnE.mean)],[],2); 
kapp_max.reverse = -min([data.V.mean ./ exp(data.lnE.mean)],[],2);
kapp_max.forward(kapp_max.forward<0) = nan;
kapp_max.reverse(kapp_max.reverse<0) = nan;

calculation_time = toc(tstart_model_balancing);


% -----------------------------------------------
% Graphics
% -----------------------------------------------

if cmb_options.show_graphics,
  display(' '); ca;
  [res, res_true, res_true_data] = cmb_graphics(prior, data, optimal, true, cmb_options, q_info, filenames.graphics_dir, kapp_max);
else
  res = []; res_true = []; res_true_data = [];
end


% -----------------------------------------------
% Save results as SBtab files and .mat files
% -----------------------------------------------

if cmb_options.save_results,
  display(' '); 
  
  cmb_save_results(network, data, bounds, optimal, filenames, cmb_options, struct('calculation_time', calculation_time, 'consistent', 1), true);
  
  if length(filenames.results_mat)
    save(filenames.results_mat,'optimal', 'calculation_time', 'gradient_V', 'init', 'cmb_options', 'V', 'kapp_max', 'preposterior', 'pp', 'filenames', 'cmb_options', 'network', 'q_info', 'prior', 'bounds', 'data', 'true', 'res', 'res_true', 'res_true_data');
  end
  
end


% -----------------------------------------------
% Result statistics
% -----------------------------------------------

if cmb_options.display + cmb_options.save_results,
  display(' '); 
  report = cmb_statistics(network,cmb_options,prior,data,optimal,calculation_time, res, res_true);
end

if cmb_options.save_results,
  if ~isempty(filenames.report_txt),
    mytable(report,0,filenames.report_txt);
  end
end

% --------------------------------------------------------------
%% Clear global variables
clearvars -global global_structure_matrices Mplus Mminus Wplus Wminus nm nr ind_M ind_Wp ind_Wm LP_info
% --------------------------------------------------------------
