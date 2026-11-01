% -------------------------------------------------------------------
% Demo script "Model balancing - test with artificial data"
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
%   network_file    (filename)        network file in SBML format(.xml) or SBtab format (.tsv)
%   constraint_file (SBtab filename): metabolite bounds    (optional; set to [] to use default constraints)
%   position_file   (SBtab filename): coordinates for network graphics (optional, otherwise set to []); needed to
%                                     define layout for SBML network files; overrides a Position table in SBtab network file
%
% Other settings:
%   model           (string) model name (can be freely chosen)
%   prb             (string) estimation scenario name (chosen by the user, should not refer to optimisation settings)
%   run             (string) estimation scenario name (chosen by the user, may refer to optimisation settings)
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
prb             = 'test';
run             = 'test_alpha_0_5';
flag_artificial = 1;
network_file    = [mb_resourcedir '/models/branch_point_model/branch_point_model.xml'];
%network_file   = [mb_resourcedir '/models/branch_point_model/branch_point_model.tsv'];
constraint_file = [mb_resourcedir '/models/branch_point_model/branch_point_model_ConcentrationConstraint.tsv'];
position_file   = [mb_resourcedir '/models/branch_point_model/branch_point_model_Position.tsv'];
c_init          = [1,1,0.1,0.5]'; 
ns              = 6;
alpha           = 0.5;
result_dir      = tempdir;

% -------------------------------------------------------------
% Small test model "three_chain_model" (network structure from SBML or SBtab file)
%
% Network structure
% X - O - O - X
%
% -------------------------------------------------------------

% model           = 'three_chain_model';
% prb             = 'test';
% run             = 'test_alpha_0_5';
% flag_artificial = 1;
% network_file    = [mb_resourcedir '/models/three_chain_model/three_chain_model.xml'];
% %network_file   = [mb_resourcedir '/models/three_chain_model/three_chain_model.tsv'];
% constraint_file = [mb_resourcedir '/models/three_chain_model/three_chain_model_ConcentrationConstraint.tsv'];
% position_file   = [mb_resourcedir '/models/three_chain_model/three_chain_model_Position.tsv'];
% c_init          = [10,1,1,0.01]'; 
% ns              = 6;
% alpha           = 0.5;
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

% model           = 'double_branch_model';
% prb             = 'test';
% run             = 'test_alpha_0_5';
% flag_artificial = 1;
% network_file    = [mb_resourcedir '/models/double_branch_model/double_branch_model.xml'];
% %network_file   = [mb_resourcedir '/models/double_branch_model/double_branch_model.tsv'];
% constraint_file = [mb_resourcedir '/models/double_branch_model/double_branch_model_ConcentrationConstraint.tsv'];
% position_file   = [mb_resourcedir '/models/double_branch_model/double_branch_model_Position.tsv'];
% c_init          = [10,1,1,1,1,1,0.01]'; 
% ns              = 6;
% alpha           = 0.5;
% result_dir      = tempdir;


% -------------------------------------------------------------
% E. coli CCM model (network structure from SBML file)
% (calculation time about 30 minutes)
% -------------------------------------------------------------

% To use this example model, please uncomment the following lines

% model           = 'e_coli_artificial';  
% prb             = 'test';
% run             = 'test_alpha_0_5';
% flag_artificial = 1;
% result_dir      = tempdir;
% network_file    = [mb_resourcedir '/models/e_coli_noor_2016/e_coli_noor_2016.tsv'];
% constraint_file = [mb_resourcedir '/models/e_coli_noor_2016/e_coli_noor_2016_ConcentrationConstraint.tsv'];
% position_file   = [mb_resourcedir '/models/e_coli_noor_2016/e_coli_noor_2016_Position.tsv'];
% c_init          = mb_e_coli_c_init(model,mb_default_options);
% alpha           = 0.5;
% ns              = 4;


% -------------------------------------------------------------
% Generate default filenames for output files
% -------------------------------------------------------------

filenames = mb_filenames(model, prb, run, result_dir, network_file,flag_artificial);


% -------------------------------------------------------------
% Set options (for other options, see 'help mb_default_options')
% -------------------------------------------------------------

mb_options                        = mb_default_options(flag_artificial);
mb_options.run                    = run;
mb_options.ns                     = ns;
mb_options.prior_variant          = 'broad_prior';
mb_options.use_kinetic_data       = 'all';
mb_options.use_artificial_noise   = 0;
mb_options.use_kinetic_data_noise = 0;
mb_options.enzyme_score_alpha     = alpha; 
mb_options.initial_values_variant = 'average_sample';
mb_options.show_graphics          = 1;
mb_options.verbose                = 0;


% -----------------------------------------------
% Prepare data structures: load SBML network model and 
% generate true kinetic model, priors, and artificial data
% -----------------------------------------------

[network, bounds, prior, q_info, data, true, kinetic_data, state_data] = mb_model_artificial_data(network_file, mb_options, c_init, position_file, constraint_file);

save_kinetic_and_state_data(network, kinetic_data, state_data, filenames.kinetic_data, filenames.state_data, 0);


% -----------------------------------------------
% Run model balancing and save results to files 
% -----------------------------------------------

model_balancing(filenames, mb_options, network, q_info, prior, bounds, data, true);
