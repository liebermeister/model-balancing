function mb_options = mb_default_options(flag_artificial)

% mb_options = mb_default_options(flag_artificial)
%
% flag_artificial (Boolean, default 0) options for model with artificial data
%                 (currently without effect; kept for compatibility)
%
% Returns matlab struct with default options for model balancing
%
% The default values are read from two files:
%  - the SBtab Config file 'model-balancing_Config.tsv' (in the main code folder 'model-balancing';
%    TableID 'ConfigureModelBalancing'): all options that do not depend on the prior table
%  - the prior table (see 'mb_prior_file'): bounds and priors in mb_options.quantities, and
%    the options metabolic_prior_c_geom_mean, metabolic_prior_c_geom_std, metabolic_prior_e_geom_mean,
%    metabolic_prior_e_geom_std, data_C_geom_std, data_E_geom_std
% To change the defaults, edit these files. The quantities can also be modified by the user in mb_options.quantities.
%
% Algorithm settings
%
%   .run                     % name of estimation scenario (freely chosen by the user)
%   .kinetic_data_set        % name of kinetic data set    (freely chosen by the user)
%   .ecm_score               % ECM score function, see 'help ecm_scores'
%   .initial_values_variant  % strategy for choosing initial values (see below and 'help mb_estimation')
%   .enzyme_score_type       % 'quadratic' (alpha=1), 'monotonic' (alpha =0), 'interpolated'
%                              (only "monotonic" guarantees that MB is convex!)
%   .enzyme_score_alpha      % stringency parameter alpha for enzyme_score_type='interpolated'
%   .beta_ln_c_over_km       % beta = 1/sigma for penalty on sum_til ln(c_i(t)/km_il) = 1/ln(geom std of c/KM)!!
%   .use_kinetic_data        % 'all', 'only_Keq', 'none'
%   .score                   % options: 'neg_log_posterior', 'log_neg_log_posterior';
%   .use_gradient            % flag: use analytic gradients in fmincon (see mb_log_posterior_gradient;
%                              only for ecm_score = 'emc4cm', otherwise numerical gradients)
%   .use_safe_optimisation   % flag: run several rounds of optimisation (until convergence)
%   .random_seed             % random seed for optimisation
%
% Display and output options
%
%   .verbose                 % flag: verbose
%   .optim_display           % option 'Display' of fmincon
%   .display                 % flag: display table with scores after each optimisation round
%   .show_graphics           % flag: show graphics
%   .plot_true_vs_data       % flag: show extra plots: true values against data values
%                              (only in the case of artificial data)
%   .save_results            % flag: save results to files
%   .save_graphics           % flag: save graphics to files
%
% Options for generating artificial data (see 'help mb_model_artificial_data')
%
%   .use_artificial_noise    % flag: generate artificial state data with noise
%   .use_kinetic_data_noise  % flag: generate artificial kinetic data with noise
%   .prior_variant           % prior variants (see mb_model_artificial_data): 'original_prior',
%                            % 'broad_prior', 'broad_prior_around_zero',
%                            % 'prior_around_true_values', 'broad_prior_around_true_values'
%
% Bounds and distributions for model variables
%   (default values from the prior file, see 'mb_prior_file', except for
%    quantities.v.max, quantities.Aforward.min, and quantities.mu0.std, which are taken from the Config file)
%
%   .use_bounds                  Boolean -> if set to 0, lower and upper bounds are ignored for all variables
%   .use_pseudo_values           Boolean -> if set to 0, no pseudo values are used
%   .use_crosscovariances        Boolean -> if set to 0, crosscovariances (between types of kinetic constants) are ignored
%   .quantities.Keq.min          new value for lower_bound(ind.Keq);
%   .quantities.Keq.max          new value for upper_bound(ind.Keq);
%   .quantities.Keq.mean_ln      new value for prior_median_ln(ind.Keq);
%   .quantities.Keq.std_ln       new value for prior_geostd_ln(ind.Keq);
%   .quantities.KV.min           new value for lower_bound(ind.KV);
%   .quantities.KV.max           new value for upper_bound(ind.KV);
%   .quantities.KV.mean_ln       new value for prior_median_ln(ind.KV);
%   .quantities.KV.std_ln        new value for prior_std_ln(ind.KV);
%   .quantities.KM.min           new value for lower_bound(ind.KM);
%   .quantities.KM.max           new value for upper_bound(ind.KM);
%   .quantities.KM.mean_ln       new value for prior_median_ln(ind.KM);
%   .quantities.KM.std_ln        new value for prior_geostd_ln(ind.KM);
%   .quantities.KA.min           new value for lower_bound(ind.KA);
%   .quantities.KA.max           new value for upper_bound(ind.KA);
%   .quantities.KA.mean_ln       new value for prior_median_ln(ind.KA);
%   .quantities.KA.std_ln        new value for prior_geostd_ln(ind.KA);
%   .quantities.KI.min           new value for lower_bound(ind.KI);
%   .quantities.KI.max           new value for upper_bound(ind.KI);
%   .quantities.KI.mean_ln       new value for prior_median_ln(ind.KI);
%   .quantities.KI.std_ln        new value for prior_geostd_ln(ind.KI);
%   .quantities.Kcatf.min        new value for lower_bound(ind.Kcatf);
%   .quantities.Kcatf.max        new value for upper_bound(ind.Kcatf);
%   .quantities.Kcatf.mean_ln    new value for prior_median_ln(ind.Kcatf);
%   .quantities.Kcatf.std_ln     new value for prior_geostd_ln(ind.Kcatf);
%   .quantities.Kcatr.min        new value for lower_bound(ind.Kcatr);
%   .quantities.Kcatr.max        new value for upper_bound(ind.Kcatr);
%   .quantities.Kcatr.mean_ln    new value for prior_median_ln(ind.Kcatr);
%   .quantities.Kcatr.std_ln     new value for prior_geostd_ln(ind.Kcatr);
%   .quantities.v.max            new value for upper_bound(v) (mM/s)
%   .quantities.c.min            new value for lower_bound(ind.c);
%   .quantities.c.max            new value for upper_bound(ind.c);
%   .quantities.e.max            new value for upper_bound(ind.u);
%   .quantities.Aforward.min     new value for lower_bound(ind.A) (kJ/mol)
%   .quantities.Aforward.max     new value for upper_bound(ind.A);
%
% Geometric standard deviations for priors
%   (also used for sampling "true" values, in the case of artificial data)
%   default values are taken from the prior file, see 'mb_prior_file'
%   .metabolic_prior_c_geom_mean
%   .metabolic_prior_e_geom_mean
%   .metabolic_prior_c_geom_std
%   .metabolic_prior_e_geom_std
%
% Geometric standard deviations for data values (describing measurement uncertainties)
%   - used in 'state_data_to_data', to complete standard deviations not given in the data files
%   - also used in 'mb_generate_artificial_data', for generating noise in artificial data
%   .data_C_geom_std     (default from prior table)
%   .data_E_geom_std     (default from prior table)
%   .data_V_geom_std     (currently not used)
%   .data_kin_geom_std
%
% Distributions for "true values" in artificial data (in 'mb_generate_artificial_data')
%   .metabolic_artificial_c_geom_std
%   .metabolic_artificial_e_geom_std
%  also used in 'mb_generate_artificial_data':
%   .quantities.mu0.std   (RT * ln 5)
%   .quantities.KV.mean_ln
%   .quantities.KV.std_ln
%   .quantities.KV.min
%   .quantities.KV.max
%   .quantities.KM.mean_ln
%   .quantities.KM.std_ln
%   .quantities.KM.min
%   .quantities.KM.max
%
% Possible initial value variants (option 'mb_options.initial_values_variant'):
%
%  'polytope center'   center of mass of some random LP solutions
%  'preposterior_mode' use preposterior mode for kinetic constants and metabolite levels,
%                      moved slightly into the interior of the allowed region
%  'flat_objective'    pick a point in the solution space by using matlab's fmincon with a uniform objective function,
%                      starting from the polytope center
%  'random'            initialise parameter vector with random values
%  'true_values'       true values: only with artificial data
%  'given_values'      use existing mb_options.init
%  'average_sample'    first run a model balancing with concentrations and fluxes averaged over all samples;
%                      and using initial value option "preposterior_mode"; use result to initialise values
%                      (if this point is not feasible, 'preposterior_mode' is used instead)
%
% Note that artificial data are only generated for ecm_score = 'emc4cm'!

eval(default('flag_artificial','0'));


% ---------------------------------------------------------
% Options from the Config file

config_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'model-balancing_Config.tsv');
mb_options  = load_config_table(config_file);


% ----------------------------------------------------------
% Read prior distributions and bounds from parameter balancing prior file

parameter_prior_file = mb_prior_file;
verbose              = 0;
parameter_prior      = parameter_balancing_prior([],parameter_prior_file,verbose);

upper_bound     = cell_string2num(parameter_prior.UpperBound);
lower_bound     = cell_string2num(parameter_prior.LowerBound);
prior_median_ln = log(cell_string2num(parameter_prior.PriorMedian));
prior_geostd_ln = log(cell_string2num(parameter_prior.PriorGeometricStd));
data_geom_std   = cell_string2num(parameter_prior.DataGeometricStd);

for it = 1:length(parameter_prior.Symbol),
  ind.(parameter_prior.Symbol{it}) = it;
end


% ----------------------------------------------------------
% Prior distributions and bounds from the prior file

for my_quantity = {'Keq','KV','KM','KA','KI','Kcatf','Kcatr'}
  q = my_quantity{1};
  mb_options.quantities.(q).min     = lower_bound(ind.(q));
  mb_options.quantities.(q).max     = upper_bound(ind.(q));
  mb_options.quantities.(q).mean_ln = prior_median_ln(ind.(q));
  mb_options.quantities.(q).std_ln  = prior_geostd_ln(ind.(q));
end

mb_options.quantities.c.min           = lower_bound(ind.c);
mb_options.quantities.c.max           = upper_bound(ind.c);
mb_options.quantities.c.mean_ln       = prior_median_ln(ind.c);
mb_options.quantities.c.std_ln        = prior_geostd_ln(ind.c);
mb_options.quantities.c.data_geom_std = data_geom_std(ind.c);

mb_options.quantities.e.max           = upper_bound(ind.u);
mb_options.quantities.e.mean_ln       = prior_median_ln(ind.u);
mb_options.quantities.e.std_ln        = prior_geostd_ln(ind.u);
mb_options.quantities.e.data_geom_std = data_geom_std(ind.u);

% (quantities.Aforward.min is taken from the Config file)
mb_options.quantities.Aforward.max    = upper_bound(ind.A);

mb_options.metabolic_prior_c_geom_mean = exp(mb_options.quantities.c.mean_ln);   % mM
mb_options.metabolic_prior_c_geom_std  = exp(mb_options.quantities.c.std_ln);    % dimensionless
mb_options.metabolic_prior_e_geom_mean = exp(mb_options.quantities.e.mean_ln);   % mM
mb_options.metabolic_prior_e_geom_std  = exp(mb_options.quantities.e.std_ln);    % dimensionless


% ---------------------------------------------------------
% Assumed distributions of measurement errors
% (= measurement noise assumed in fitting procedure,
%  = noise level used when generating artificial data with noise)
% data_V_geom_std and data_kin_geom_std are taken from the Config file

mb_options.data_C_geom_std = mb_options.quantities.c.data_geom_std;
mb_options.data_E_geom_std = mb_options.quantities.e.data_geom_std;
