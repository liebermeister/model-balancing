% -------------------------------------------------------------------
% Convex parameter estimation demo script
% 
% This demo script is a simplified version of the script cmb_example_artificial_data.m
% -------------------------------------------------------------------

model      = 'small_model'; 
run        = 'test';
result_dir = '/tmp/';


% -------------------------------------------------------------
% Define model and initial concentration guess
% -------------------------------------------------------------

network_sbml = [cmb_resourcedir '/models/small_model/small_model.xml'];
c_init       = [10,1,1,1,1,1,0.01]'; 


% -------------------------------------------------------------
% Set default options
% -------------------------------------------------------------

cmb_options                        = cmb_default_options;
cmb_options.run                    = run;
cmb_options.prior_variant          = 'broad_prior';
cmb_options.use_artificial_noise   = 1;
cmb_options.use_kinetic_data       = 'all';
cmb_options.initial_values_variant = 'average_sample'; 
cmb_options.show_graphics          = 1;


% -------------------------------------------------------------
% Set run-specific options
% -------------------------------------------------------------

cmb_options.ns                     = 6; 
cmb_options.use_kinetic_data       = 'all';
cmb_options.use_artificial_noise   = 0;
cmb_options.use_kinetic_data_noise = 0;
cmb_options.initial_values_variant = 'average_sample';


% -----------------------------------------------
% Prepare data structures: load SBML network model 
% and generate "true" kinetic model, priors, and artificial data
% -----------------------------------------------

[filenames, network, q_info, prior, bounds, data, true] = cmb_model_artificial_data(model, network_sbml, cmb_options, c_init, result_dir);


% -----------------------------------------------
% Run convex parameter estimation
% -----------------------------------------------

convex_model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, true);
