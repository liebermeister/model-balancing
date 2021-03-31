% -------------------------------------------------------------
% Model balancing - Simple test model
%
% o model and data are read from an SBtab file
% o model balancing is performed
% -------------------------------------------------------------

clear;

% -------------------------------------------------------------
% Define model
% -------------------------------------------------------------

% -------------------------------------------------------------
% Small test model "branch_point_model" (network structure from SBML or SBtab file)
%
% Network structure
%
% X
%  \   
%    O - X
%  /  
% X
% 
% Use pregenerated artificial data from folder ~/resources/models/branch_point_model/data/
% -------------------------------------------------------------

model       = 'branch_point_model'; 
run         = 'with_kinetic_data_balanced';
result_dir  = tempdir;
filenames   = cmb_filenames(model, run, result_dir);
sbtab_files = {[cmb_basedir '/../resources/models/branch_point_model/data/artificial_network_true.tsv'], ...
               [cmb_basedir '/../resources/models/branch_point_model/data/artificial_data.tsv']};
file_options.match_data_by = 'ModelElementId';
file_options.columns_mean  = {'S1_Mean','S2_Mean','S3_Mean','S4_Mean','S5_Mean','S6_Mean'};
file_options.columns_std   = {'S1_Std','S2_Std','S3_Std','S4_Std','S5_Std','S6_Std'};


% -------------------------------------------------------------
% Small test model "three_chain_model" (network structure from SBML or SBtab file)
%
% Network structure
%
% X - O - O - X
%
% Use pregenerated artificial data from folder ~/resources/models/three_chain_model/data/
% -------------------------------------------------------------

% model       = 'three_chain_model'; 
% run         = 'with_kinetic_data_balanced';
% result_dir  = tempdir;
% filenames   = cmb_filenames(model, run, result_dir);
% sbtab_files = {[cmb_basedir '/../resources/models/three_chain_model/data/artificial_network_true.tsv'], ...
%                [cmb_basedir '/../resources/models/three_chain_model/data/artificial_data.tsv']};
% file_options.match_data_by = 'ModelElementId';
% file_options.columns_mean  = {'S1_Mean','S2_Mean','S3_Mean','S4_Mean','S5_Mean','S6_Mean'};
% file_options.columns_std   = {'S1_Std','S2_Std','S3_Std','S4_Std','S5_Std','S6_Std'};


% -------------------------------------------------------------
% Small test model "double_branch_model" (network structure from SBML or SBtab file)
%
% Network structure
%
% X - O       X
%      \     /
%       O - O 
%      /     \
% X - O       X
%
% Use pregenerated artificial data from folder ~/resources/models/double_branch_model/data/
% -------------------------------------------------------------

% model       = 'double_branch_model'; 
% run         = 'with_kinetic_data_balanced';
% result_dir  = tempdir;
% filenames   = cmb_filenames(model, run, result_dir);
% sbtab_files = {[cmb_basedir '/../resources/models/double_branch_model/data/artificial_network_true.tsv'], ...
%                [cmb_basedir '/../resources/models/double_branch_model/data/artificial_data.tsv']};
% file_options.match_data_by = 'ModelElementId';
% file_options.columns_mean  = {'S1_Mean','S2_Mean','S3_Mean','S4_Mean','S5_Mean','S6_Mean'};
% file_options.columns_std   = {'S1_Std','S2_Std','S3_Std','S4_Std','S5_Std','S6_Std'};


% -------------------------------------------------
% Set options
% -------------------------------------------------

cmb_options                          = cmb_default_options;
cmb_options.run                      = run;
cmb_options.prior_variant            = 'broad_prior';
cmb_options.initial_values_variant   = 'average_sample'; 
cmb_options.metabolic_prior_geom_std = 20;
cmb_options.save_results             = 0;

switch run, 
  
  case 'with_kinetic_data_original',    
    cmb_options.use_kinetic_data       = 'all';
    cmb_options.source_of_kinetic_data = 'original';
  
  case 'with_Keq_data_original',
    cmb_options.use_kinetic_data       = 'only_Keq_data';
    cmb_options.source_of_kinetic_data = 'original';
  
  case 'no_kinetic_data_original',
    cmb_options.use_kinetic_data       = 'none';
    cmb_options.source_of_kinetic_data = 'original';
  
  case 'with_kinetic_data_balanced',    
    cmb_options.use_kinetic_data       = 'all';
    cmb_options.source_of_kinetic_data = 'balanced';
  
  case 'with_Keq_data_balanced',
    cmb_options.use_kinetic_data       = 'only_Keq_data';
    cmb_options.source_of_kinetic_data = 'balanced';
  
  case 'no_kinetic_data_balanced',
    cmb_options.use_kinetic_data       = 'none';
    cmb_options.source_of_kinetic_data = 'balanced';

end


% -----------------------------------------------
% Load model
% -----------------------------------------------

[network, q_info, data, c_init, kinetic_data, state_data] = cmb_load_and_convert_model_and_data(cmb_options, sbtab_files, file_options);


% -----------------------------------------------
% Generate kinetic model and priors
% -----------------------------------------------

[bounds, prior, init] = cmb_model_and_data(model, network, data, q_info, c_init, cmb_options);

% Test: show imported data 

if 0,
  kinetic_data_print(kinetic_data,network);
  [c_mean, c_std ] = lognormal_log2normal(data.X.mean, data.X.std);
  [c_mean, c_std ]
  [data.E.mean, data.E.std]
  [data.V.mean, data.V.std]
end


% --------------------------------------------------------------
% Display data (metabolite levels and fluxes, glucose condition)
% --------------------------------------------------------------

netgraph_concentrations(network, network.external, data.V.mean(:,1),1); % exp(data.X.mean(:,1))


% --------------------------------------------------------------
% Save input data 
% --------------------------------------------------------------

if cmb_options.save_results,
  save_kinetic_and_state_data(network, kinetic_data, state_data, filenames.kinetic_data, filenames.state_data);
end


% -----------------------------------------------
% Run convex parameter estimation
% -----------------------------------------------

model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, [], init);
