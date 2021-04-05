function cmb_options = cmb_default_options()

% cmb_options = cmb_default_options()
%
% Returns matlab struct with default options for model balancing
% Fields and default values are shown below
%
% Algorithm settings
% 
%   .run                     = 'my_run';            % name of estimation scenario (freely chosen by the user)
%   .kinetic_data_set        = 'original';          % name of kinetic data set    (freely chosen by the user)
%   .ecm_score               = 'emc4cm';            % ECM score function, see 'help ecm_scores'
%   .initial_values_variant  = 'average_sample';    % strategy for choosing initial values (see 'help cmb_estimation')
%                                                     possible choices: 'average_sample', 'preposterior_mode', 
%                                                     'random', 'true_values', 'given_values'
%   .enzyme_score_type       = 'interpolated';      % 'quadratic' (alpha=1), 'monotonic' (alpha =0), 'interpolated'
%                                                     (only "monotonic" guarantees that MB is convex!)
%   .enzyme_score_alpha      = 0.5;                 % interpolation parameter for enzyme_score_type='interpolated'
%   .parameterisation        = 'Keq_KV_KM_KA_KI';   % options: 'Keq_KV_KM_KA_KI', 'KV_KM_KA_KI';
%   .use_kinetic_data        = 'all';               % 'all', 'only_Keq_data', 'none'
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
%   .prior_variant           = 'original_prior';                original prior
%                              'broad_prior'                    broad prior around original prior means
%                              'broad_prior_around_zero'        broad prior around 0 values
%                              'prior_around_true_values'       prior around true values (only for artificial data)
%                              'broad_prior_around_true_values' broad prior around true values (only for artificial data)
%
% Bounds and distributions for model variables
%   (Options to override the default values set in the prior table file)
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
%   .quantities.v.max            new value for upper_bound(v), default 100 mM/s
%   .quantities.c.min            new value for lower_bound(ind.c);
%   .quantities.c.max            new value for upper_bound(ind.c);
%   .quantities.e.max            new value for upper_bound(ind.u);
%   .quantities.Aforward.min     new value for lower_bound(ind.A), default 0.0001 kJ/Mol
%   .quantities.Aforward.max     new value for upper_bound(ind.A);
% 
% Distributions from which "true values" are drawn
%
%   .quantities.mu0.std            default RT * log(5);
%                                  
%   .metabolic_prior_c_geom_mean   default 1 mM
%   .metabolic_prior_e_geom_mean   default 0.001 mM
%   .metabolic_prior_c_geom_std    default 10
%   .metabolic_prior_e_geom_std    default 10
%                                  
%   .data_kin_geom_std             default 1.5
%   .data_V_geom_std               default 1.2 
%   .data_C_geom_std               default 1.2 
%   .data_E_geom_std               default 1.2 
% 
% Distributions for noise in artificial data
% 
%   .metabolic_artificial_c_geom_std   default 1.5;
%   .metabolic_artificial_e_geom_std   default 1.5;
% 


% ---------------------------------------------------------
% Algorithm settings

cmb_options.run                    = 'my_run';
cmb_options.ecm_score              = 'emc4cm'; 
% note that artificial data are only generated for ecm_score = 'emc4cm'!


% Variants: 'preposterior_mode', 'random', 'true_values', 'given_values', 'average_sample'
cmb_options.initial_values_variant = 'average_sample';
cmb_options.enzyme_score_type      = 'interpolated'; % 'monotonic'; % 'quadratic' (only "monotonic" guarantees that MB is convex!)
cmb_options.enzyme_score_alpha     = 0.5; % 'quadratic' (only "monotonic" guarantees that MB is convex!)

cmb_options.use_kinetic_data       = 'all';
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
%cmb_options.quantities.KV.std_ln    = prior_geostd_ln(ind.KV);
cmb_options.quantities.KV.std_ln     = log(5);

cmb_options.quantities.KM.min       = lower_bound(ind.KM);
cmb_options.quantities.KM.max       = upper_bound(ind.KM);
cmb_options.quantities.KM.mean_ln   = prior_median_ln(ind.KM);
% cmb_options.quantities.KM.std_ln    = prior_geostd_ln(ind.KM);
cmb_options.quantities.KM.std_ln    = log(5);
 
cmb_options.quantities.KA.min       = lower_bound(ind.KA);
cmb_options.quantities.KA.max       = upper_bound(ind.A);
cmb_options.quantities.KA.mean_ln   = prior_median_ln(ind.KA);
cmb_options.quantities.KA.std_ln    = prior_geostd_ln(ind.KA);

cmb_options.quantities.KI.min       = lower_bound(ind.KI);
cmb_options.quantities.KI.max       = upper_bound(ind.KI);
cmb_options.quantities.KI.mean_ln   = prior_median_ln(ind.KI);
cmb_options.quantities.KI.std_ln    = prior_geostd_ln(ind.KI);

cmb_options.quantities.v.max        = 100;                % mM/s

cmb_options.quantities.c.min        = lower_bound(ind.c);
cmb_options.quantities.c.max        = upper_bound(ind.c);

cmb_options.quantities.e.max        = upper_bound(ind.u);

cmb_options.quantities.Aforward.min = 0.0001; % kJ/Mol
%cmb_options.quantities.Aforward.min = 0.1; % kJ/Mol
cmb_options.quantities.Aforward.max = upper_bound(ind.A);

cmb_options.metabolic_prior_c_geom_mean   = 1;     % mM
cmb_options.metabolic_prior_e_geom_mean   = 0.001; % mM
cmb_options.metabolic_prior_c_geom_std      = 10;    % 
cmb_options.metabolic_prior_e_geom_std      = 10;    % 


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
cmb_options.data_kin_geom_std = 1.5;

cmb_options.data_V_geom_std   = 1.2; 
cmb_options.data_C_geom_std   = 1.2; 
cmb_options.data_E_geom_std   = 1.2; 

% ---------------------------------------------------------
% distributions for artificial data

% distributions from which "true values" are drawn

cmb_options.quantities.mu0.std = RT * log(5);

% distributions for noise in artificial data

cmb_options.metabolic_artificial_c_geom_std = 1.5;
cmb_options.metabolic_artificial_e_geom_std = 1.5;

