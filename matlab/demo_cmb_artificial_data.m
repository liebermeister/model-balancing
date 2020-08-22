% -------------------------------------------------------------------
% Convex parameter estimation demo script
% -------------------------------------------------------------------


% -------------------------------------------------------------
% Define model and initial concentration guess
% -------------------------------------------------------------

% Small test model

model        = 'small_model'; 
run          = 'test';
result_dir   = tempdir;
network_sbml = [cmb_resourcedir '/models/small_model/small_model.xml'];
c_init       = [10,1,1,1,1,1,0.01]'; 
ns           = 6;

% E. coli CCM model - to balance this model, uncomment the following lines
% (calculation time about 30 minutes)

% model        = 'e_coli_artificial';  
% run          = 'no_noise_with_noisefree_kinetic_data';
% result_dir   = tempdir;
% network_sbml = [cmb_resourcedir '/data/data-organisms/escherichia_coli/network/ecoli_ccm_ProteinComposition_Network.xml'];
% c_init       = e_coli_c_init(model,cmb_default_options);
% ns           = 4;


% -------------------------------------------------------------
% Set options
% -------------------------------------------------------------

cmb_options                        = cmb_default_options;
cmb_options.run                    = run;
cmb_options.ns                     = ns;
cmb_options.prior_variant          = 'broad_prior';
cmb_options.use_kinetic_data       = 'all';
cmb_options.use_artificial_noise   = 0;
cmb_options.use_kinetic_data_noise = 0;
cmb_options.initial_values_variant = 'average_sample';
cmb_options.show_graphics          = 1;
cmb_options.verbose                = run;


% -----------------------------------------------
% Prepare data structures: load SBML network model and generate "true" 
% kinetic model, priors, and artificial data
% -----------------------------------------------

filenames = cmb_filenames(model, run, result_dir, network_sbml);

[network, bounds, prior, q_info, data, true, kinetic_data, state_data] = cmb_model_artificial_data(model, network_sbml, cmb_options, c_init);


% -----------------------------------------------
% Run convex parameter estimation
% -----------------------------------------------

convex_model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, true);
