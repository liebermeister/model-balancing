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
% ~/resources/models/double_branch_model/data, used in demo_mb_experimental_data.m
%
% -------------------------------------------------------------

% -------------------------------------------------------------
% MOST OF THE CODE BELOW STEMS FROM demo_mb_artificial_data.m
% -------------------------------------------------------------

model           = 'double_branch_model';
run             = 'test';
result_dir      = tempdir;
network_file    = [mb_resourcedir '/models/double_branch_model/double_branch_model.xml'];
position_file   = [mb_resourcedir '/models/double_branch_model/double_branch_model_Position.tsv'];
constraint_file = [mb_resourcedir '/models/double_branch_model/double_branch_model_ConcentrationConstraint.tsv'];
c_init          = [10,1,1,1,1,1,0.01]'; 
ns              = 6;

filenames = mb_filenames(model, run, result_dir, network_file);

filenames.network_true = [mb_resourcedir '/models/double_branch_model/data/artificial_network_true.tsv'];
filenames.kinetic_data = [mb_resourcedir '/models/double_branch_model/data/artificial_kinetic_data.tsv'];
filenames.state_data = [mb_resourcedir '/models/double_branch_model/data/artificial_state_data.tsv'];


% -------------------------------------------------------------
% Set options
% -------------------------------------------------------------

mb_options                        = mb_default_options;
mb_options.run                    = run;
mb_options.ns                     = ns;
mb_options.prior_variant          = 'broad_prior';
mb_options.use_kinetic_data       = 'all';
mb_options.use_artificial_noise   = 0;
mb_options.use_kinetic_data_noise = 0;
mb_options.initial_values_variant = 'average_sample';
mb_options.show_graphics          = 1;
mb_options.verbose                = 0;


% -----------------------------------------------
% Prepare data structures: load SBML network model and 
% generate "true" kinetic model, priors, and artificial data
% -----------------------------------------------

[network, bounds, prior, q_info, data, true_model, kinetic_data, state_data] = mb_model_artificial_data(network_file, mb_options, c_init, position_file, constraint_file);

save_kinetic_and_state_data(network, kinetic_data, state_data, filenames.kinetic_data, filenames.state_data, 0);

network.kinetics = true_model.kinetics; 
network_to_sbtab(network,struct('filename', filenames.network_true));
