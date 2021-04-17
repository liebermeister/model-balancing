function cmb_options = cmb_default_options(flag_artificial)

% cmb_options = cmb_default_options(flag_artificial)
%
% flag_artificial (Boolean, default 0) options for model with artificial data
% 
% Returns matlab struct with default options for model balancing
% Fields and default values are shown below
% Note that this function reads information about priors from the prior table (see 'cmb_prior_file')
% and stores them in cmb_options.quantities; there they can be modified by the user
%
% Algorithm settings
% 
%   .run                     = 'my_run';            % name of estimation scenario (freely chosen by the user)
%   .kinetic_data_set        = 'original';          % name of kinetic data set    (freely chosen by the user)
%   .ecm_score               = 'emc4cm';            % ECM score function, see 'help ecm_scores'
%   .initial_values_variant  = 'average_sample';    % strategy for choosing initial values (see below and 'help cmb_estimation')
%                                                     possible choices: 'average_sample', 'preposterior_mode', 
%                                                     'random', 'true_values', 'given_values'
%   .enzyme_score_type       = 'interpolated';      % 'quadratic' (alpha=1), 'monotonic' (alpha =0), 'interpolated'
%                                                     (only "monotonic" guarantees that MB is convex!)
%   .enzyme_score_alpha      = 0.5;                 % stringency parameter alpha for enzyme_score_type='interpolated'
%   .beta_ln_c_over_km       = 0                    % beta = 1/sigma for penalty on sum_til ln(c_i(t)/km_il)
%   .parameterisation        = 'Keq_KV_KM_KA_KI';   % options: 'Keq_KV_KM_KA_KI', 'KV_KM_KA_KI';
%   .use_kinetic_data        = 'all';               % 'all', 'only_Keq', 'none'
%   .score                   = 'neg_log_posterior'; % options: 'neg_log_posterior', 'log_neg_log_posterior';
%   .use_gradient            = 0;                   % flag: set opt.SpecifyObjectiveGradient = true in fmincon optimisation
%   .use_safe_optimisation   = 1;                   % flag: run several rounds of optimisation (until convergence)
%   .random_seed             = 1;                   % random seed for optimisation
% 
% Display and output options
%
%   .verbose                 = 0;                   % flag: verbose
%   .optim_display           = 'off';               % flag: show graphical output during optimisation 
%   .display                 = 1;                   % flag: display table with scores after each optimisation round
%   .show_graphics           = 1;                   % flag: show graphics
%   .plot_true_vs_data       = 0;                   % flag: show extra plots: true values against data values 
%                                                     (only in the case of artificial data)
%   .save_results            = 1;                   % flag: save results to files
%   .save_graphics           = 1;                   % flag: save graphics to files
%
% Options for generating artificial data (see 'help cmb_model_artificial_data')
%
%   .use_artificial_noise    = 0;                   % flag: generate artificial state data with noise
%   .use_kinetic_data_noise  = 1;                   % flag: generate artificial kinetic data with noise
%   .prior_variant           = 'original_prior';    % prior variants (see cmb_model_artificial_data)
%                                                   % also 'broad_prior', 'broad_prior_around_zero', 
%                                                   % 'prior_around_true_values', 'broad_prior_around_true_values'
%
% Bounds and distributions for model variables
%   Default values are taken from the prior file, see 'cmb_prior_file'
%
%   .quantities.Keq.max          new value for upper_bound(ind.Keq);
%   .quantities.Keq.mean_ln      new value for prior_median_ln(ind.Keq);
%   .quantities.Keq.std_ln       new value for prior_geostd_ln(ind.Keq);
%   .quantities.KV.min           new value for lower_bound(ind.KV);
%   .quantities.KV.max           new value for upper_bound(ind.KV);
%   .quantities.KV.mean_ln       new value for prior_median_ln(ind.KV);
%   .quantities.KV.std_ln        new value for prior_std_ln(ind.KV), default log(5);
%   .quantities.KM.min           new value for lower_bound(ind.KM);
%   .quantities.KM.max           new value for upper_bound(ind.KM);
%   .quantities.KM.mean_ln       new value for prior_median_ln(ind.KM);
%   .quantities.KM.std_ln        new value for prior_geostd_ln(ind.KM), default log(5 mM)
%   .quantities.KA.min           new value for lower_bound(ind.KA);
%   .quantities.KA.max           new value for upper_bound(ind.A);
%   .quantities.KA.mean_ln       new value for prior_median_ln(ind.KA);
%   .quantities.KA.std_ln        new value for prior_geostd_ln(ind.KA);
%   .quantities.KI.min           new value for lower_bound(ind.KI);
%   .quantities.KI.max           new value for upper_bound(ind.KI);
%   .quantities.KI.mean_ln       new value for prior_median_ln(ind.KI);
%   .quantities.KI.std_ln        new value for prior_geostd_ln(ind.KI);
%   .quantities.Kcatf.min           new value for lower_bound(ind.Kcatf);
%   .quantities.Kcatf.max           new value for upper_bound(ind.Kcatf);
%   .quantities.Kcatf.mean_ln       new value for prior_median_ln(ind.Kcatf);
%   .quantities.Kcatf.std_ln        new value for prior_geostd_ln(ind.Kcatf);
%   .quantities.Kcatr.min           new value for lower_bound(ind.Kcatr);
%   .quantities.Kcatr.max           new value for upper_bound(ind.Kcatr);
%   .quantities.Kcatr.mean_ln       new value for prior_median_ln(ind.Kcatr);
%   .quantities.Kcatr.std_ln        new value for prior_geostd_ln(ind.Kcatr);
%   .quantities.v.max            new value for upper_bound(v), default 100 mM/s
%   .quantities.c.min            new value for lower_bound(ind.c);
%   .quantities.c.max            new value for upper_bound(ind.c);
%   .quantities.e.max            new value for upper_bound(ind.u);
%   .quantities.Aforward.min     new value for lower_bound(ind.A), default 0.0001 kJ/Mol
%   .quantities.Aforward.max     new value for upper_bound(ind.A);
% 
% Geometric standard deviations for priors 
%   (also used for sampling "true" values, in the case of artificial data)
%   default values are taken from the prior file, see 'cmb_prior_file'
%   .metabolic_prior_c_geom_mean   
%   .metabolic_prior_e_geom_mean   
%   .metabolic_prior_c_geom_std    
%   .metabolic_prior_e_geom_std    
%                           
% Geometric standard deviations for data values (describing measurement uncertainties)
%   - used in 'state_data_to_data', to complete standard deviations not given in the data files
%   - also used in 'cmb_generate_artificial_data', for generating noise in artificial data
%   .data_C_geom_std     default from prior table
%   .data_E_geom_std     default from prior table
%   .data_V_geom_std     default 1.2 % currently not used
%   .data_kin_geom_std   default 1.5 
% 
% Distributions for "true values" inartificial data (in 'cmb_generate_artificial_data')
%   .metabolic_artificial_c_geom_std   default 1.5;
%   .metabolic_artificial_e_geom_std   default 1.5;
%  also used in 'cmb_generate_artificial_data':
%   .quantities.mu0.std                default RT * log(5);
%   .quantities.KV.mean_ln
%   .quantities.KV.std_ln
%   .quantities.KV.min
%   .quantities.KV.max
%   .quantities.KM.mean_ln
%   .quantities.KM.std_ln
%   .quantities.KM.min
%   .quantities.KM.max
%
% Possible initial value variants (option 'cmb_options.initial_values_variant'):
% 
%  'polytope center'   (default) center of mass of some random LP solutions
%  'preposterior_mode' use preposterior mode for kinetic constants and metabolite levels
%  'flat_objective'    pick a point in the solution space by using matlab's fmincon with a uniform objective function
%  'random'            initialise parameter vector with random values
%  'true_values'       true values: only with artificial data
%  'given_values'      use existing cmb_options.init
%  'average_sample'    first run a model balancing with concentrations and fluxes averaged over all samples; 
%                      and using initial value option "preposterior_mode"; use result to initialise values

eval(default('flag_artificial','0'));

  
% ---------------------------------------------------------
% Algorithm settings

cmb_options.run                    = 'my_run';
cmb_options.ecm_score              = 'emc4cm'; 
% note that artificial data are only generated for ecm_score = 'emc4cm'!


% Variants: 'preposterior_mode', 'random', 'true_values', 'given_values', 'average_sample'
cmb_options.initial_values_variant = 'polytope center';
cmb_options.enzyme_score_type      = 'interpolated'; % 'monotonic'; % 'quadratic' (only "monotonic" guarantees that MB is convex!)
cmb_options.enzyme_score_alpha     = 0.5; % 'quadratic' (only "monotonic" guarantees that MB is convex!)

cmb_options.use_kinetic_data       = 'all';
cmb_options.beta_ln_c_over_km      = 0;
cmb_options.parameterisation       = 'Keq_KV_KM_KA_KI'; % 'KV_KM_KA_KI';


cmb_options.use_artificial_noise   = 0; % use noise in artificial state data 
cmb_options.use_kinetic_data_noise = 1; % use noise in artificial kinetic data 

% used to modify prior for tests with artificial data (see cmb_model_artificial_data.m)
cmb_options.prior_variant          = 'original_prior'; 
cmb_options.use_safe_optimisation  = 1;
cmb_options.use_gradient           = 0;

% also 'log_neg_log_posterior', see cmb_objective.m:
cmb_options.score                  = 'neg_log_posterior'; 
cmb_options.kinetic_data_set = 'original'; % balanced

cmb_options.random_seed       = 1;
cmb_options.verbose           = 0;
cmb_options.optim_display     = 'off';
cmb_options.display           = 1;
cmb_options.show_graphics     = 1;
cmb_options.plot_true_vs_data = 0;
cmb_options.save_results      = 1;
cmb_options.save_graphics     = 1;
  
% ----------------------------------------------------------
% Read prior distributions and bounds from parameter balancing prior file

parameter_prior_file = cmb_prior_file; 
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
% Prior distributions and bounds (mostly from parameter balancing prior file)

cmb_options.quantities.Keq.max      = upper_bound(ind.Keq);
cmb_options.quantities.Keq.mean_ln  = prior_median_ln(ind.Keq);
cmb_options.quantities.Keq.std_ln   = prior_geostd_ln(ind.Keq);

cmb_options.quantities.KV.min       = lower_bound(ind.KV);
cmb_options.quantities.KV.max       = upper_bound(ind.KV);
cmb_options.quantities.KV.mean_ln   = prior_median_ln(ind.KV);
cmb_options.quantities.KV.std_ln    = prior_geostd_ln(ind.KV);

cmb_options.quantities.KM.min       = lower_bound(ind.KM);
cmb_options.quantities.KM.max       = upper_bound(ind.KM);
cmb_options.quantities.KM.mean_ln   = prior_median_ln(ind.KM);
cmb_options.quantities.KM.std_ln    = prior_geostd_ln(ind.KM);
 
cmb_options.quantities.KA.min       = lower_bound(ind.KA);
cmb_options.quantities.KA.max       = upper_bound(ind.A);
cmb_options.quantities.KA.mean_ln   = prior_median_ln(ind.KA);
cmb_options.quantities.KA.std_ln    = prior_geostd_ln(ind.KA);

cmb_options.quantities.KI.min       = lower_bound(ind.KI);
cmb_options.quantities.KI.max       = upper_bound(ind.KI);
cmb_options.quantities.KI.mean_ln   = prior_median_ln(ind.KI);
cmb_options.quantities.KI.std_ln    = prior_geostd_ln(ind.KI);

cmb_options.quantities.Kcatf.min     = lower_bound(ind.Kcatf);
cmb_options.quantities.Kcatf.max     = upper_bound(ind.Kcatf);
cmb_options.quantities.Kcatf.mean_ln = prior_median_ln(ind.Kcatf);
cmb_options.quantities.Kcatf.std_ln  = prior_geostd_ln(ind.Kcatf);

cmb_options.quantities.Kcatr.min     = lower_bound(ind.Kcatr);
cmb_options.quantities.Kcatr.max     = upper_bound(ind.Kcatr);
cmb_options.quantities.Kcatr.mean_ln = prior_median_ln(ind.Kcatr);
cmb_options.quantities.Kcatr.std_ln  = prior_geostd_ln(ind.Kcatr);

cmb_options.quantities.v.max        = 100; % mM/s

cmb_options.quantities.c.min        = lower_bound(ind.c);
cmb_options.quantities.c.max        = upper_bound(ind.c);
cmb_options.quantities.c.mean_ln    = prior_median_ln(ind.c);
cmb_options.quantities.c.std_ln     = prior_geostd_ln(ind.c);
cmb_options.quantities.c.data_geom_std = data_geom_std(ind.c);

cmb_options.quantities.e.max        = upper_bound(ind.u);
cmb_options.quantities.e.mean_ln    = prior_median_ln(ind.u);
cmb_options.quantities.e.std_ln     = prior_geostd_ln(ind.u);
cmb_options.quantities.e.data_geom_std = data_geom_std(ind.u);

cmb_options.quantities.Aforward.min = 0.0001; % kJ/Mol
%cmb_options.quantities.Aforward.min = 0.1;    % kJ/Mol
%cmb_options.quantities.Aforward.min = lower_bound(ind.A);
cmb_options.quantities.Aforward.max = upper_bound(ind.A);

cmb_options.metabolic_prior_c_geom_mean   = exp(cmb_options.quantities.c.mean_ln);   % mM
cmb_options.metabolic_prior_c_geom_std    = exp(cmb_options.quantities.c.std_ln);    % dimensionless 
cmb_options.metabolic_prior_e_geom_mean   = exp(cmb_options.quantities.e.mean_ln);   % mM
cmb_options.metabolic_prior_e_geom_std    = exp(cmb_options.quantities.e.std_ln);    % dimensionless

% ----------------------------------------------------
% Manual choice

% cmb_options.quantities.Keq.max      = 100000;
% cmb_options.quantities.Keq.mean_ln  = log(1);
% cmb_options.quantities.Keq.std_ln   = log(10);
% cmb_options.quantities.KV.min       = 0.00001;
% cmb_options.quantities.KV.max       = 10000;
% cmb_options.quantities.KV.max       = 10000;
% cmb_options.quantities.KM.min       = 0.0001;
% cmb_options.quantities.KA.min       = 0.001;
% cmb_options.quantities.KI.min       = 0.001;
% cmb_options.quantities.KV.mean_ln   = log(20);
% cmb_options.quantities.KV.std_ln     = log(5);
% cmb_options.quantities.KM.max       = 100;
% cmb_options.quantities.KM.mean_ln   = log(1);
% cmb_options.quantities.KM.std_ln    = log(5);
% cmb_options.quantities.KA.max       = 100;
% cmb_options.quantities.KA.mean_ln   = log(1);
% cmb_options.quantities.KA.std_ln    = log(10);
% cmb_options.quantities.KI.max       = 100;
% cmb_options.quantities.KI.mean_ln   = log(1);
% cmb_options.quantities.KI.std_ln    = log(10);
% cmb_options.quantities.c.min        = 0.000001;
% cmb_options.quantities.c.max        = 100;
% cmb_options.quantities.v.max        = 100;
% cmb_options.quantities.e.max        = 10;
% cmb_options.quantities.Aforward.max = 50;        % kJ/Mol
% cmb_options.quantities.Aforward.min = 0.0000001; % kJ/Mol


% ---------------------------------------------------------
% assumed distributions of measurement errors

% = measurement noise assumed in fitting procedure
% = noise level used when generating artificial data with noise

% relatively high noise level, assuming deviation between in vitro and in vivo data 

  cmb_options.data_C_geom_std   = cmb_options.quantities.c.data_geom_std;
  cmb_options.data_E_geom_std   = cmb_options.quantities.e.data_geom_std; 
  cmb_options.data_V_geom_std   = 1.2;
  cmb_options.data_kin_geom_std = 1.5;

% ---------------------------------------------------------
% distributions for artificial data

if flag_artificial,

  %% distributions from which "true values" are drawn
  
  cmb_options.quantities.mu0.std = RT * log(5);
  
  %% distributions for noise in artificial data
  
  cmb_options.metabolic_artificial_c_geom_std = 1.5;
  cmb_options.metabolic_artificial_e_geom_std = 1.5;

end
