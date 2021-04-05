% -------------------------------------------------------------------
% Demo script "Model balancing - test with articifial data"
%
% What this script does:
% o read model from SBML or SBtab file
% o generate artificial data (kinetic constants and state data)
% o run model balancing to try and recover the true model and state data
% -------------------------------------------------------------------


% -------------------------------------------------------------
% Define model and options
% -------------------------------------------------------------

% Settings: 
%
% Files that define a model: 
%   network_file    (SBML format(.xml) or SBtab format (.tsv))
%   position_file   (SBtab): graphics coordinates (optional; set to [] to use default positions)
%   constraint_file (SBtab): metabolite bounds    (optional; set to [] to use default constraints)
%
% Other settings:
%   model           (string) model name (can be freely chosen)
%   run             (string) estimation scenario name (can be freely chosen)
%   c_init          (column vector) metabolite concentration profile (optional; set to [] to use default vector of 1's) 
%   ns              number of metabolic states
%   result_dir      (directory name) directory for output files (default = matlab's tempdir) 


% -------------------------------------------------------------
% Small test model "branch_point_model" (network structure from SBML or SBtab file)
%
% Network structure
% X
%  \   
%    O - X
%  /  
% X
% -------------------------------------------------------------

model           = 'branch_point_model';
run             = 'test';
network_file    = [cmb_resourcedir '/models/branch_point_model/branch_point_model.xml'];
%network_file   = [cmb_resourcedir '/models/branch_point_model/branch_point_model.tsv'];
position_file   = [cmb_resourcedir '/models/branch_point_model/branch_point_model_Position.tsv'];
constraint_file = [cmb_resourcedir '/models/branch_point_model/branch_point_model_ConcentrationConstraint.tsv'];
c_init          = [1,1,0.1,0.5]'; 
ns              = 6;
result_dir      = tempdir;


% -------------------------------------------------------------
% Small test model "three_chain_model" (network structure from SBML or SBtab file)
%
% Network structure
% X - O - O - X
%
% -------------------------------------------------------------

% model           = 'three_chain_model';
% run             = 'test';
% network_file   = [cmb_resourcedir '/models/three_chain_model/three_chain_model.xml'];
% %network_file    = [cmb_resourcedir '/models/three_chain_model/three_chain_model.tsv'];
% position_file   = [cmb_resourcedir '/models/three_chain_model/three_chain_model_Position.tsv'];
% constraint_file = [cmb_resourcedir '/models/three_chain_model/three_chain_model_ConcentrationConstraint.tsv'];
% c_init          = [10,1,1,1,1,1,0.01]'; 
% ns              = 6;
% result_dir      = tempdir;


% -------------------------------------------------------------
% Small test model "double_branch_model" (network structure from SBML or SBtab file)
%
% Network structure
% X - O       X
%      \     /
%       O - O 
%      /     \
% X - O       X
% -------------------------------------------------------------

%  model        = 'double_branch_model';
%  run          = 'test';
%  network_file   = [cmb_resourcedir '/models/double_branch_model/double_branch_model.xml'];
%  %network_file    = [cmb_resourcedir '/models/double_branch_model/double_branch_model.tsv'];
%  position_file   = [cmb_resourcedir '/models/double_branch_model/double_branch_model_Position.tsv'];
%  constraint_file = [cmb_resourcedir '/models/double_branch_model/double_branch_model_ConcentrationConstraint.tsv'];
%  c_init       = [10,1,1,1,1,1,0.01]'; 
%  ns           = 6;
%  result_dir   = tempdir;


% -------------------------------------------------------------
% E. coli CCM model (network structure from SBML file)
% (calculation time about 30 minutes)
% -------------------------------------------------------------

% To use this example model, please uncomment the following lines

% model         = 'e_coli_artificial';  
% run           = 'test';
% result_dir    = tempdir;
% network_file  = [cmb_resourcedir '/models/e_coli_noor_2016/ecoli_ccm_ProteinComposition_Network.xml'];
% position_file = [];
% c_init        = cmb_e_coli_c_init(model,cmb_default_options);
% ns            = 4;


% -------------------------------------------------------------
% Generate default filenames for output files
% -------------------------------------------------------------

filenames = cmb_filenames(model, run, result_dir, network_file);


% -------------------------------------------------------------
% Set options (for other options, see 'help cmb_default_options')
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
% generate true kinetic model, priors, and artificial data
% -----------------------------------------------

[network, bounds, prior, q_info, data, true, kinetic_data, state_data] = cmb_model_artificial_data(network_file, cmb_options, c_init, position_file, constraint_file);

save_kinetic_and_state_data(network, kinetic_data, state_data, filenames.kinetic_data, filenames.state_data, 0);


% -----------------------------------------------
% Run model balancing and save results to files 
% -----------------------------------------------

model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, true);
