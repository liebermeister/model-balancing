% -------------------------------------------------------------
% Model "double_branch_model" (network structure from SBML or SBtab file)
%
% Network structure
% X - O       X
%      \     /
%       O - O 
%      /     \
% X - O       X
%
% This script generates the artificial model and data files in 
% ~/resources/models/double_branch_model/data, used in demo_cmb_experimental_data.m
%
% -------------------------------------------------------------

% -------------------------------------------------------------
% MOST OF THE CODE BELOW STEMS FROM demo_cmb_artificial_data.m
% -------------------------------------------------------------

model           = 'double_branch_model';
run             = 'test';
result_dir      = tempdir;
network_file    = [cmb_resourcedir '/models/double_branch_model/double_branch_model.xml'];
position_file   = [cmb_resourcedir '/models/double_branch_model/double_branch_model_Position.tsv'];
constraint_file = [cmb_resourcedir '/models/double_branch_model/double_branch_model_ConcentrationConstraint.tsv'];
c_init          = [10,1,1,1,1,1,0.01]'; 
ns              = 6;

filenames = cmb_filenames(model, run, result_dir, network_file);

filenames.network_true = '/home/wolfram/projekte/model_balancing/github/model-balancing/resources/models/double_branch_model/data/artificial_network_true.tsv';
filenames.kinetic_data = '/home/wolfram/projekte/model_balancing/github/model-balancing/resources/models/double_branch_model/data/artificial_kinetic_data.tsv';
filenames.state_data = '/home/wolfram/projekte/model_balancing/github/model-balancing/resources/models/double_branch_model/data/artificial_state_data.tsv';


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
cmb_options.verbose                = 0;


% -----------------------------------------------
% Prepare data structures: load SBML network model and 
% generate "true" kinetic model, priors, and artificial data
% -----------------------------------------------

[network, bounds, prior, q_info, data, true, kinetic_data, state_data] = cmb_model_artificial_data(network_file, cmb_options, c_init, position_file, constraint_file);

save_kinetic_and_state_data(network, kinetic_data, state_data, filenames.kinetic_data, filenames.state_data, 0);

network.kinetics = true.kinetics; 
network_to_sbtab(network,struct('filename', filenames.network_true));
