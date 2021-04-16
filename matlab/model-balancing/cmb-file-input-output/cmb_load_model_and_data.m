function [network, kinetic_data, state_data, conc_min, conc_max] = cmb_load_model_and_data(model_file, kinetic_data_file, state_data_files, replace_ids_in_network, match_data_by, constraint_sbtab_file)

% [network, kinetic_data, state_data, conc_min, conc_max] = cmb_load_model_and_data(model_file, kinetic_data_file, state_data_files, replace_ids_in_network, match_data_by, constraint_sbtab_file)
%
% Load model and (kinetic and state) data from (SBML and SBtab) files into matlab data structures
%   1. load a model structure (from sbml or sbtab file) (possibly including default metabolite levels c_init)
%   2. load kinetic data and metabolic data (from sbtab files) and map them to the model via KEGG IDs
%
% Input:
%   model_file:           model file (SBML or SBtab format)
%   kinetic_data_file:    one (filename) or several (list of filenames) kinetic data files (SBtab format)
%                         instead of a filename, one may directly rovide a kinetic data structure 
%                         from which data will be selected
%   metabolite_data_file: metabolite data (filename)
%   enzyme_data_file:     enzyme data (filename)
%   flux_data_file:       flux data (filename)
%
%   replace_ids_in_network = [];
%   match_data_by          = 'KeggId'; % ModelElementId
%
% Output:
%   network:        model data structure (as in MNT toolbox)
%   kinetic_data:   kinetic data data structure (as in MNT toolbox, mnt_kinetic_data)
%   state_data:     metabolite, enzyme, and flux data

eval(default('kinetic_data_file', '[]', 'metabolite_data_file', '[]', 'enzyme_data_file', '[]', 'flux_data_file', '[]','replace_ids_in_network','[]','match_data_by','''KeggId'''));
  
options.use_sbml_ids         = 0;
options.use_kegg_ids         = 1;
options.parameter_prior_file = cmb_prior_file;
  
display(sprintf('Importing model and data from file %s',model_file));

network          = network_import_model(model_file, struct('load_quantity_table',0));
network.kinetics = set_kinetics(network,'cs');

if ~isfield(network,'metabolite_KEGGID'), network.metabolite_KEGGID = {}; end
if ~isfield(network,'reaction_KEGGID'),   network.reaction_KEGGID = {};   end

if length(kinetic_data_file),
  parameter_prior = parameter_balancing_prior([],options.parameter_prior_file,0); 
  [model_quantities, basic_quantities, data_quantities, pseudo_quantities] = parameter_balancing_quantities(parameter_prior, network, options);
  if isempty(find(cellfun('length',network.metabolite_KEGGID))), warning('No KEGG species annotations found in model');  end
  if isempty(find(cellfun('length',network.reaction_KEGGID))),   warning('No KEGG species annotations found in model');  end
  opt = struct('use_sbml_ids', options.use_sbml_ids, 'use_kegg_ids', options.use_kegg_ids, 'verbose', 0);
  kinetic_data = kinetic_data_load(data_quantities, parameter_prior, network, kinetic_data_file, opt);
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
  sdf.metabolite.replace_ids_in_network = replace_ids_in_network;
  sdf.enzyme.replace_ids_in_network     = replace_ids_in_network;
  sdf.flux.replace_ids_in_network       = replace_ids_in_network;
  sdf.metabolite.match_data_by          = match_data_by;
  sdf.enzyme.match_data_by              = match_data_by;
  sdf.flux.match_data_by                = match_data_by;

  state_data.metabolite_data = load_state_data_for_network(network, sdf.metabolite.file, sdf.metabolite.type, sdf.metabolite, sdf.metabolite.table);
  state_data.enzyme_data     = load_state_data_for_network(network, sdf.enzyme.file, sdf.enzyme.type, sdf.enzyme, sdf.enzyme.table);
  state_data.flux_data       = load_state_data_for_network(network, sdf.flux.file, sdf.flux.type, sdf.flux, sdf.flux.table);
  if isfield(sdf,'delta_g'), 
    state_data.delta_g_data  = load_state_data_for_network(network, sdf.delta_g.file, sdf.delta_g.type, sdf.delta_g, sdf.delta_g.table);
  end
end

% project fluxes or complete missing fluxes by projections
% NOTE THAT PROJECTION REQUIRES ASSUMES THAT ALL FLUXES MUST BE IN FORWARD DIRECTION

if find(isnan(state_data.flux_data.Mean)),
  display('- (cmb_load_model_and_data.m): projecting fluxes to complete missing values');
  flag_project_fluxes = 'missing'; % 'all', 'missing', 'none'
  [state_data.flux_data.Mean, state_data.flux_data.Std] = cmb_project_fluxes(state_data.flux_data.Mean,state_data.flux_data.Std,network,flag_project_fluxes);
end


% concentration bounds
conc_min = [];
conc_max = [];
if length(constraint_sbtab_file),
  constraint_sbtab = sbtab_document_load_from_one(constraint_sbtab_file);
  conc_min  = sbtab_table_get_column(constraint_sbtab.tables.ConcentrationConstraint,'Min',1);
  conc_max  = sbtab_table_get_column(constraint_sbtab.tables.ConcentrationConstraint,'Max',1);
end
