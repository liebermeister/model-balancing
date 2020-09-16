% -------------------------------------------------------------
% Model balancing - Simple test model
%
% o model and data are read from an SBtab file
% o model balancing is performed
% -------------------------------------------------------------

% -------------------------------------------------------------
% Simple test model
% 
% Network structure
% -\   /
%    -
% _/   \
% -------------------------------------------------------------

clear;

model      = 'double_branch_model'; 
run        = 'with_kinetic_data_balanced';
result_dir = tempdir;
filenames  = cmb_filenames(model, run, result_dir);
sbtab_dir  = [cmb_basedir '/../resources/models/double_branch_model/sbtab'];
sbtab_file = [sbtab_dir '/double_branch_model.tsv'];

%% Alternative: use data from individual SBtab files
%% file_network    = [dir '/multi/double_branch_model_NETWORK.tsv'];
%% file_kinetic    = [dir '/multi/double_branch_model_KINETIC.tsv'];
%% file_metabolite = [dir '/multi/double_branch_model_METABOLITE.tsv'];
%% file_flux       = [dir '/multi/double_branch_model_FLUX.tsv'];
%% file_enzyme     = [dir '/multi/double_branch_model_ENZYME.tsv'];
%% sbtab_file      = {file_network, file_kinetic, file_metabolite, file_flux, file_enzyme};


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

[network, q_info, data, c_init, kinetic_data, state_data] = cmb_read_model_and_data(cmb_options, sbtab_file);


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

convex_model_balancing(filenames, cmb_options, network, q_info, prior, bounds, data, [], init);
