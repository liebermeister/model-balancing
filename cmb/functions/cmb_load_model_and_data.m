function [network, kinetic_data, state_data] = cmb_load_model_and_data(model_file, kinetic_data_file, state_data_files, replace_ids_in_network)

% [network, kinetic_data, state_data] = cmb_load_model_and_data(model_file, kinetic_data_file, metabolite_data_file, enzyme_data_file, flux_data_file)
%
% Load model and (kinetic and state) data from (SBML and SBtab) files
%   1. load a model structure (from sbml or sbtab file) (possibly including default metabolite levels c_init)
%   2. load kinetic data and metabolic data (from sbtab files) and map them to the model via KEGG IDs
%
% Input:
%   model_file:           model file (SBML or SBtab format)
%   kinetic_data_file:    one (filename) or several (list of filenames) kinetic data files (SBtab format)
%                         instead of a filename, one may directly rovide a kinetic data structure from which data will be selected
%   metabolite_data_file: metabolite data (filename)
%   enzyme_data_file:     enzyme data (filename)
%   flux_data_file:       flux data (filename)
%
% Output:
%   network:        model data structure (as in MNT toolbox)
%   kinetic_data:   kinetic data data structure (as in MNT toolbox, mnt_kinetic_data)
%   state_data:     metabolite, enzyme, and flux data


eval(default('kinetic_data_file', '[]', 'metabolite_data_file', '[]', 'enzyme_data_file', '[]', 'flux_data_file', '[]','replace_ids_in_network','[]'));
  
options.use_sbml_ids         = 0;
options.use_kegg_ids         = 1;
options.parameter_prior_file = cmb_prior_file;
  
network = network_import_model(model_file);

if length(kinetic_data_file),
  parameter_prior = parameter_balancing_prior([],options.parameter_prior_file,0); 
  [model_quantities, basic_quantities, data_quantities, pseudo_quantities] = parameter_balancing_quantities(parameter_prior, network, options);
  if isempty(find(cellfun('length',network.metabolite_KEGGID))), warning('No KEGG species annotations found in model');  end
  if isempty(find(cellfun('length',network.reaction_KEGGID))),   warning('No KEGG species annotations found in model');  end
  opt = struct('use_sbml_ids', options.use_sbml_ids, 'use_kegg_ids', options.use_kegg_ids, 'verbose', 0);
  kinetic_data = data_integration_load_kinetic_data(data_quantities, parameter_prior, network, kinetic_data_file, opt);
else
  kinetic_data = [];
end

metabolite_data = [];
enzyme_data     = [];
flux_data       = [];

if length(state_data_files),
  sdf = state_data_files;
  if isfield(state_data_files,'samples'),
    state_data.samples = state_data_files.samples;
  end
  state_data.metabolite_data = load_network_state_data(network, sdf.metabolite.file, sdf.metabolite.type, sdf.metabolite.columns_mean, sdf.metabolite.columns_std, replace_ids_in_network);
  state_data.enzyme_data     = load_network_state_data(network, sdf.enzyme.file, sdf.enzyme.type, sdf.enzyme.columns_mean, sdf.enzyme.columns_std, replace_ids_in_network);
  state_data.flux_data       = load_network_state_data(network, sdf.flux.file, sdf.flux.type, sdf.flux.columns_mean, sdf.flux.columns_std, replace_ids_in_network);
end
