function [network, q_info, data, c_init, kinetic_data, state_data, conc_min, conc_max] = cmb_load_and_convert_model_and_data(cmb_options, file, file_options)

% [network, q_info, data, c_init, kinetic_data, state_data, conc_min, conc_max] = cmb_load_and_convert_model_and_data(cmb_options, file, file_options)
%
% Load model and data from SBtab files for model balancing
% 
% Wrapper function around:
%  - sbtab_load_network_model_and_data
%  - state_data_to_data
%
% This function does not accept SBML files as inputs (to be implemented)
%
% Function arguments:
%  'cmb_options':      data structure with model balancing options; 
%                        for list of options used in this function, see state_data_to_data
%  'file':             either a single SBtab filename or a cell array of such filenames (files to be combined)
%  'file_options':     struct with (optional) fields
%    .match_data_by         'ModelElementId' or 'KeggId';
%    .metabolite_table_id   (default 'MetaboliteConcentration')
%    .flux_table_id         (default 'Flux')
%    .enzyme_table_id       (default 'EnzymeConcentration')
%    .columns_mean:         cell array of column names for mean data values
%                           (same columns names in metabolite, flux, and enzyme table!)
%    .columns_std:          same, for std dev columns
%    .columns_std:          same, for std dev columns
%    .constraint_sbtab_file (optional) SBtab file with metabolite concentration constraints 
%                           (table 'ConcentrationConstraint')
%
% Outputs
% conc_min, conc_max        vectors of concentration bounds (or empty, if no constraint_sbtab_file is given) 
  
eval(default('file_options','struct'));

file_options = join_struct(struct('constraint_sbtab_file',[]), file_options);

% --------------------------------------------------------------
% Load data files -> network, kinetic_data, state_data
% --------------------------------------------------------------

result = sbtab_load_network_model_and_data(file, file_options);

network          = result.network;
network.kinetics = set_kinetics(network,'cs');
kinetic_data     = result.kinetic_data;
state_data       = result.state_data;

% --------------------------------------------------------------
% Initial state for optimisation (using parameters from network.kinetics and concentrations = 1)
% --------------------------------------------------------------

nm = size(network.N,1);
c_init = ones(nm,1);


% --------------------------------------------------------------
% Convert metabolic data to struct 'data'
% --------------------------------------------------------------

[data, q_info] = state_data_to_data(kinetic_data, state_data, cmb_options, network);


% --------------------------------------------------------------
% concentration bound vectors
% --------------------------------------------------------------

conc_min = [];
conc_max = [];

if length(file_options.constraint_sbtab_file),
  constraint_sbtab = sbtab_document_load_from_one(file_options.constraint_sbtab_file);
  conc_min  = sbtab_table_get_column(constraint_sbtab.tables.ConcentrationConstraint,'Min',1);
  conc_max  = sbtab_table_get_column(constraint_sbtab.tables.ConcentrationConstraint,'Max',1);
end
