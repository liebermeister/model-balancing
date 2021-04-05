% -------------------------------------------------------------
% Model balancing - Simple test models
%
% What this scrip does:
% o read model and data from SBtab file
% o run model balancing
% o save results to output files
%
% -------------------------------------------------------------

clear;

% -------------------------------------------------------------
% Define model
% -------------------------------------------------------------

% Settings:
%
% Files that define model and data: 
%   network_file       (SBtab filename) SBtab file with tables 'Reaction', 'Compound', 
%                                          'Position' (optional), 'Parameter' (optional)
%   kinetic_data_file  (SBtab filename) SBtab file with table  'ParameterData' (optional; if kinetic_data_file == [],
%                                          kinetic data are read from table 'Parameter' in file network_file)
%   state_data_file    (SBtab filename) SBtab file with tables 'MetabolicFluxData', 
%                                          'MetaboliteConcentrationData', 'EnzymeConcentrationData'
%
% Other settings:
%   model                        (string) model name (freely chosen by the user)
%   run                          (string) estimation scenario name (freely chosen by the user)
%   file_options.match_data_by   'ModelElementId' (default) or 'KeggId'
%   file_options.columns_mean    (cell array of strings) column names for mean data values in metabolite, flux, and enzyme table
%   file_options.columns_std     (cell array of strings) column names for std deviations   in metabolite, flux, and enzyme table
%   file_options.constraint_file (SBtab filename): file with metabolite bounds (optional; set to [] to use default constraints)
%   result_dir                   (directory name) directory for output files (default = matlab's tempdir) 


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
% Use model and artificial data from folder ~/resources/models/branch_point_model/data/
% -------------------------------------------------------------

model                        = 'branch_point_model'; 
run                          = 'with_kinetic_data';
network_file                = [cmb_resourcedir '/models/branch_point_model/data/artificial_network_true.tsv'];
kinetic_data_file           = [cmb_resourcedir '/models/branch_point_model/data/artificial_kinetic_data.tsv'];
state_data_file             = [cmb_resourcedir '/models/branch_point_model/data/artificial_state_data.tsv'];
file_options.constraint_file = [cmb_resourcedir '/models/branch_point_model/branch_point_model_ConcentrationConstraint.tsv'];
file_options.match_data_by   = 'ModelElementId';
file_options.columns_mean    = {'S1_Mean','S2_Mean','S3_Mean','S4_Mean','S5_Mean','S6_Mean'};
file_options.columns_std     = {'S1_Std','S2_Std','S3_Std','S4_Std','S5_Std','S6_Std'};
result_dir                   = tempdir; % replace by your desired output directory

filenames = cmb_filenames(model, run, result_dir);


% -------------------------------------------------------------
% Small test model "three_chain_model" (network structure from SBML or SBtab file)
%
% Network structure
%
% X - O - O - X
%
% Use model and artificial data from folder ~/resources/models/three_chain_model/data/
% -------------------------------------------------------------

% model                      = 'three_chain_model'; 
% run                        = 'with_kinetic_data';
% network_file              = [cmb_resourcedir '/models/three_chain_model/data/artificial_network_true.tsv'];
% kinetic_data_file         = [cmb_resourcedir '/models/three_chain_model/data/artificial_kinetic_data.tsv'];
% state_data_file           = [cmb_resourcedir '/models/three_chain_model/data/artificial_state_data.tsv'];
% file_options.constraint_file = [cmb_resourcedir '/models/branch_point_model/three_chain_model_ConcentrationConstraint.tsv'];
% file_options.match_data_by = 'ModelElementId';
% file_options.columns_mean  = {'S1_Mean','S2_Mean','S3_Mean','S4_Mean','S5_Mean','S6_Mean'};
% file_options.columns_std   = {'S1_Std','S2_Std','S3_Std','S4_Std','S5_Std','S6_Std'};
% result_dir                 = tempdir; % replace by your desired output directory
%
% filenames   = cmb_filenames(model, run, result_dir);


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
% Use model and artificial data from folder ~/resources/models/double_branch_model/data/
% -------------------------------------------------------------

% model                      = 'double_branch_model'; 
% run                        = 'with_kinetic_data';
% network_file              = [cmb_resourcedir '/models/double_branch_model/data/artificial_network_true.tsv'];
% kinetic_data_file         = [cmb_resourcedir '/models/double_branch_model/data/artificial_kinetic_data.tsv'];
% state_data_file           = [cmb_resourcedir '/models/double_branch_model/data/artificial_state_data.tsv'];
% file_options.constraint_file = [cmb_resourcedir '/models/branch_point_model/double_branch_model_ConcentrationConstraint.tsv'];
% file_options.match_data_by = 'ModelElementId';
% file_options.columns_mean  = {'S1_Mean','S2_Mean','S3_Mean','S4_Mean','S5_Mean','S6_Mean'};
% file_options.columns_std   = {'S1_Std','S2_Std','S3_Std','S4_Std','S5_Std','S6_Std'};
% result_dir                 = tempdir; % replace by your desired output directory
%
% filenames   = cmb_filenames(model, run, result_dir);


% -------------------------------------------------------------
% E. coli CCM model (network structure from SBML file)
% (calculation time: a few minutes)
% -------------------------------------------------------------

% % To use this example model, please uncomment the following lines
% 
% model                     = 'e_coli_noor_2016';  
% run                       = 'test';
% network_file              = [cmb_resourcedir '/models/e_coli_noor_2016/e_coli_noor_2016.tsv'];
% kinetic_data_file         = [cmb_resourcedir '/models/e_coli_noor_2016/e_coli_ccm_kinetic_data.tsv'];
% state_data_file           = [cmb_resourcedir '/models/e_coli_noor_2016/e_coli_ccm_state_data.tsv'];
% %% Concentration constraint table is part of the model file:
% file_options.constraint_file = [cmb_resourcedir '/models/e_coli_noor_2016/ecoli_ccm_ProteinComposition_Network.xml'];
% file_options.match_data_by = 'ModelElementId';
% file_options.columns_mean  = {'S_1_Mean'};
% file_options.columns_std   = {};
% result_dir                 = tempdir; % replace by your desired output directory
% 
% filenames   = cmb_filenames(model, run, result_dir);


% -------------------------------------------------
% Set options (for other options, see 'help cmb_default_options')
% -------------------------------------------------

cmb_options                          = cmb_default_options;
cmb_options.run                      = run;
cmb_options.prior_variant            = 'broad_prior';
cmb_options.initial_values_variant   = 'average_sample'; 
cmb_options.enzyme_score_alpha       = 0.5;
cmb_options.metabolic_prior_geom_std = 10; % set width of metabolite prior
cmb_options.save_results             = 0;  % set to 1 to save result files 
cmb_options.save_graphics            = 0;  % set to 1 to save graphics files 


% -------------------------------------------------
% Set specific options for different estimation scenarios
% (please modify for your own estimation scenarios; 'run' can be freely defined)
% -------------------------------------------------

switch run, 
  
  case 'with_kinetic_data',
    cmb_options.use_kinetic_data = 'all';
  
  case 'with_Keq_data',
    cmb_options.use_kinetic_data = 'only_Keq_data';
  
  case 'no_kinetic_data',
    cmb_options.use_kinetic_data = 'none';

end


% -----------------------------------------------
% Load model and data
% -----------------------------------------------

[network, q_info, data, c_init, kinetic_data, state_data, conc_min, conc_max] = cmb_load_and_convert_model_and_data(cmb_options, {network_file, kinetic_data_file, state_data_file}, file_options);

% Save model and data to SBtab files (in result directory)

if cmb_options.save_results,
  save_kinetic_and_state_data(network, kinetic_data, state_data, filenames.kinetic_data, filenames.state_data);
end

% Please uncomment to show data

if 0,
  kinetic_data_print(kinetic_data,network);
  [c_mean, c_std ] = lognormal_log2normal(data.X.mean, data.X.std);
  [c_mean, c_std ]
  [data.E.mean, data.E.std]
  [data.V.mean, data.V.std]
end

% Please uncomment to plot network and data (mean log concentrations and fluxes for sample #1 )

if 0,
  netgraph_concentrations(network, data.X.mean(:,1), data.V.mean(:,1),1);
end


% -----------------------------------------------
% Generate data structures for bounds and priors
% -----------------------------------------------

[bounds, prior, init] = cmb_model_and_data(model, network, data, q_info, c_init, cmb_options, conc_min, conc_max);


% --------------------------------------------------------------
% Save model and data to JSON files (for python code) (in result directory)
% --------------------------------------------------------------

if cmb_options.save_results,
  cvxpy_problem = cvxpy_problem_data_structure(network, q_info, prior, data, []);
  save_json_file(cvxpy_problem, filenames.model_and_data_json);
  display(sprintf('Writing JSON file for cvxpy: %s', filenames.model_and_data_json));
end


% -----------------------------------------------
% Run model balancing and save output data 
% (file locations are specified in struct 'filenames')
% -----------------------------------------------

model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, [], init);
