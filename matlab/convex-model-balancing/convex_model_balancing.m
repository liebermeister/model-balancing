function [optimal, calculation_time, gradient_V, init, cmb_options, V, kapp_max, preposterior, pp] = convex_model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, true, init)

% [optimal, calculation_time, gradient_V, init, cmb_options, V, kapp_max, preposterior, pp] = convex_model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, true, init)
%
% Convex model balancing
%
% Input variables:
%   filenames    struct, (model- and run-specific) filenames for input and output files
%   cmb_options  struct, options 
%   network      struct, metabolic network (format as in Metabolic Network Toolbox)
%   q_info       struct, the dependencies between model variables
%   prior        struct, priors in the optimality problem (details see cmb_make_prior)
%   bounds       struct, bounds in the optimality problem (details see cmb_make_bounds)
%   data         struct, data used in the optimality problem
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

tic
    
eval(default('true','[]'));

ns = size(data.X.mean,2);
cmb_options.ns = ns;

% -----------------------------------------------
% Initialise some variables

[nr,nm] = network_numbers(network);

ns = size(data.X.mean,2);
nq = length(prior.q.mean);

% Prepare the preposterior distributions (individual terms)
preposterior = cmb_prepare_posterior(prior, data, cmb_options, q_info);

if sum(isnan(data.V.mean)),
  error('Some flux values are unknown');
end

V = data.V.mean;

pp = cmb_make_pp(network);

if isempty(true),
  if strcmp(cmb_options.initial_values_variant,'true_values'),
    warning('No true values given; using method "preposterior_mode" for choosing initial values');
    cmb_options.initial_values_variant = 'preposterior_mode';
  end
end


% -----------------------------------------------
% Initial values

switch cmb_options.initial_values_variant,
  case 'preposterior_mode';
    cmb_options.init = []; % DEFAULT: use preposterior mode (of q alone) as starting point
  case 'random',
    display('Using random initial values; please note that these values may violate the constraints');
    cmb_options.init.q = 0.1*randn(nq,1);
    cmb_options.init.X = 0.1*randn(size(data.X.mean));
  case 'true_values',
    cmb_options.init.q = true.q;
    cmb_options.init.X = true.X;
  case 'given_values';
    if ~isfield(cmb_options,'init'), 
      error('missing options field "cmb_options.init", required by "cmb_options.initial_values_variant=given_values"'); 
    end
    %% existing cmb_options.init is used 
  case 'average_sample';
    display(sprintf('\n-------------------\nFinding initial point by running model balancing on average state data'));

    my_cmb_options                        = cmb_options;
    my_cmb_options.ns                     = 1;
    my_cmb_options.initial_values_variant = 'preposterior_mode'; 
    my_cmb_options.show_graphics          = 0;
    my_cmb_options.verbose                = 0;
    my_cmb_options.display                = 0;
    my_q_info        = cmb_define_parameterisation(network, my_cmb_options); 
    [~, my_prior, my_bounds, my_data] = cmb_generate_artificial_data(network, my_cmb_options, my_q_info,ones(size(network.metabolites)));
    my_data.V.mean = nanmean(data.V.mean,2);
    my_data.X.mean = nanmean(data.X.mean,2);
    my_data.E.mean = exp(nanmean(log(data.E.mean),2));
    %% taking the average over standard deviations is not ideal .. maybe change this?
    my_data.V.std  = nanmean(data.V.std,2); 
    my_data.X.std  = nanmean(data.X.std,2);
    my_data.E.std  = exp(nanmean(log(data.E.std),2));
    my_optimal = convex_model_balancing(filenames, my_cmb_options, network, my_q_info, my_prior, my_bounds, my_data);
    cmb_options.init.q = my_optimal.q;
    cmb_options.init.X = repmat(my_optimal.X,1,cmb_options.ns);
    display(sprintf('Initial point found\n-------------------'));

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

if length(cmb_options.init),
  cmb_options.init.q(isnan(cmb_options.init.q)) = 0;
  cmb_options.init.X(isnan(cmb_options.init.X)) = 0;
end

display(' '); display('Running optimisation ..');

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

kapp_max.forward =  max([data.V.mean ./ data.E.mean],[],2); 
kapp_max.reverse = -min([data.V.mean ./ data.E.mean],[],2);
kapp_max.forward(kapp_max.forward<0) = nan;
kapp_max.reverse(kapp_max.reverse<0) = nan;

calculation_time = toc;


% -----------------------------------------------
% Graphics
% -----------------------------------------------

if cmb_options.show_graphics,
  display(' '); ca;
  cmb_graphics(prior, data, optimal, true, cmb_options, q_info, filenames.graphics_dir, kapp_max)
end


% -----------------------------------------------
% Save results as SBtab files and .mat files
% -----------------------------------------------

if cmb_options.save_results,
  display(' '); 
  
  cmb_save_results(network, data, optimal, filenames, cmb_options, struct('calculation_time',calculation_time), true);

  save(filenames.result_file,'optimal', 'calculation_time', 'gradient_V', 'init', 'cmb_options', 'V', 'kapp_max', 'preposterior', 'pp', 'filenames', 'cmb_options', 'network', 'q_info', 'prior', 'bounds', 'data', 'true');

end


% -----------------------------------------------
% Result statistics
% -----------------------------------------------

if cmb_options.display + cmb_options.save_results,
  display(' '); 
  report = cmb_statistics(network,cmb_options,prior,data,optimal,calculation_time);
end

if cmb_options.save_results,
  mytable(report,0,filenames.report);
end

