function cmb_options = cmb_default_options()

% cmb_options = cmb_default_options()
%
% Default options for convex model balancing (default values shown below)
%
% Bounds and distributions for model variables:
%
% cmb_options.quantities.Keq.max      = upper_bound(ind.Keq);
% cmb_options.quantities.Keq.mean_ln  = prior_median_ln(ind.Keq);
% cmb_options.quantities.Keq.std_ln   = prior_geostd_ln(ind.Keq);
% cmb_options.quantities.KV.min       = lower_bound(ind.KV);
% cmb_options.quantities.KV.max       = upper_bound(ind.KV);
% cmb_options.quantities.KV.mean_ln   = prior_median_ln(ind.KV);
% cmb_options.quantities.KV.std_ln    = log(5);
% cmb_options.quantities.KM.min       = lower_bound(ind.KM);
% cmb_options.quantities.KM.max       = upper_bound(ind.KM);
% cmb_options.quantities.KM.mean_ln   = prior_median_ln(ind.KM);
% cmb_options.quantities.KM.std_ln    = prior_geostd_ln(ind.KM);
% cmb_options.quantities.KM.std_ln    = log(5);
% cmb_options.quantities.KA.min       = lower_bound(ind.KA);
% cmb_options.quantities.KA.max       = upper_bound(ind.A);
% cmb_options.quantities.KA.mean_ln   = prior_median_ln(ind.KA);
% cmb_options.quantities.KA.std_ln    = prior_geostd_ln(ind.KA);
% cmb_options.quantities.KI.min       = lower_bound(ind.KI);
% cmb_options.quantities.KI.max       = upper_bound(ind.KI);
% cmb_options.quantities.KI.mean_ln   = prior_median_ln(ind.KI);
% cmb_options.quantities.KI.std_ln    = prior_geostd_ln(ind.KI);
% cmb_options.quantities.v.max        = 100;                % mM/s
% cmb_options.quantities.c.min        = lower_bound(ind.c);
% cmb_options.quantities.c.max        = upper_bound(ind.c);
% cmb_options.quantities.e.max        = upper_bound(ind.u);
% cmb_options.quantities.Aforward.min = 0.0001; % kJ/Mol
% cmb_options.quantities.Aforward.max = upper_bound(ind.A);
% 
% Distributions from which "true values" are drawn
% cmb_options.quantities.mu0.std = RT * log(5);
%
% cmb_options.metabolic_prior_c_geom_mean   = 1;     % mM
% cmb_options.metabolic_prior_e_geom_mean   = 0.001; % mM
% cmb_options.metabolic_prior_c_geom_std      = 10;    % 
% cmb_options.metabolic_prior_e_geom_std      = 10;    % 
% 
% cmb_options.data_kin_geom_std = 1.5;
% cmb_options.data_V_geom_std   = 1.2; 
% cmb_options.data_C_geom_std   = 1.2; 
% cmb_options.data_E_geom_std   = 1.2; 
% 
% Distributions for noise in artificial data
% 
% cmb_options.metabolic_artificial_c_geom_std = 1.5;
% cmb_options.metabolic_artificial_e_geom_std = 1.5;
% 
% Other settings
% 
% cmb_options.run                     = 'my_run';
% cmb_options.ecm_score               = 'emc4cm'; 
% cmb_options.optim_display           = 'off';
% cmb_options.initial_values_variant  = 'average_sample';
%                            Variants : 'preposterior_mode', 'random', 'true_values', 'given_values', 'average_sample'
% cmb_options.enzyme_likelihood_type  = 'interpolated'; % 'quadratic', 'interpolated' (only "monotonic" guarantees that MB is convex!)
% cmb_options.enzyme_likelihood_alpha = 0.5; % interpolation parameter for enzyme_likelihood_type='interpolated'
% 
% cmb_options.use_artificial_noise    = 0; % use noise in artificial data 
% cmb_options.use_kinetic_data_noise  = 1;
% cmb_options.use_kinetic_data        = 'none';
% cmb_options.parameterisation        = 'Keq_KV_KM_KA_KI'; % 'KV_KM_KA_KI';
% cmb_options.prior_variant           = 'original_prior'; 
% cmb_options.use_safe_optimisation   = 1;
% cmb_options.use_gradient            = 0;
% cmb_options.show_graphics           = 1;
% cmb_options.score                   = 'neg_log_posterior'; % 'log_neg_log_posterior';
% cmb_options.source_of_kinetic_data  = 'original'; % balanced
% cmb_options.random_seed             = 1;
% cmb_options.plot_true_vs_data       = 0;
% cmb_options.verbose                 = 0;
% cmb_options.display                 = 1;
% cmb_options.save_results            = 1;

  
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


% ---------------------------------------------------------
% Other settings

 cmb_options.run                   = 'my_run';
cmb_options.ecm_score              = 'emc4cm'; 
% note that artificial data are only generated for ecm_score = 'emc4cm'!

cmb_options.optim_display          = 'off';

% Variants: 'preposterior_mode', 'random', 'true_values', 'given_values', 'average_sample'
cmb_options.initial_values_variant = 'average_sample';
cmb_options.enzyme_likelihood_type = 'interpolated'; % 'monotonic'; % 'quadratic' (only "monotonic" guarantees that MB is convex!)
cmb_options.enzyme_likelihood_alpha= 0.5; % 'quadratic' (only "monotonic" guarantees that MB is convex!)

cmb_options.use_artificial_noise   = 0; % use noise in artificial data 
cmb_options.use_kinetic_data_noise = 1;
cmb_options.use_kinetic_data       = 'none';
cmb_options.parameterisation       = 'Keq_KV_KM_KA_KI'; % 'KV_KM_KA_KI';

% used to modify prior for tests with artificial data (see cmb_model_artificial_data.m)
cmb_options.prior_variant          = 'original_prior'; 
cmb_options.use_safe_optimisation  = 1;
cmb_options.use_gradient           = 0;
cmb_options.show_graphics          = 1;

% also 'log_neg_log_posterior', see cmb_objective.m:
cmb_options.score                  = 'neg_log_posterior'; 
cmb_options.source_of_kinetic_data = 'original'; % balanced

cmb_options.random_seed       = 1;
cmb_options.plot_true_vs_data = 0;
cmb_options.verbose           = 0;
cmb_options.display           = 1;
cmb_options.save_results      = 1;
