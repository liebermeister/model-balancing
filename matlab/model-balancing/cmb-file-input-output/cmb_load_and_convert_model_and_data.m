function [network, q_info, data, c_init, kinetic_data, state_data] = cmb_load_and_convert_model_and_data(cmb_options, file, file_options)

% [network, q_info, data, c_init, kinetic_data, state_data] = cmb_load_and_convert_model_and_data(cmb_options, file, file_options)
%
% Wrapper function around:
%  - sbtab_load_network_model_and_data
%  - state_data_to_data
  
% 'cmb_options':  see cmb_default_options.m
% 'file':         either a single sbtab filename or a cell array of such filenames (files to be combined)
% 'file_options': struct with (optional) fields
%      file_options.match_data_by       'ModelElementId' or 'KeggId';
%      file_options.metabolite_table_id (default 'MetaboliteConcentration')
%      file_options.flux_table_id       (default 'Flux')
%      file_options.enzyme_table_id     (default 'EnzymeConcentration')
%      file_options.columns_mean:       cell array of column names for mean data values
%                                       (same columns names in metabolite, flux, and enzyme table!)
%      file_options.columns_std:        same, for std dev columns

eval(default('file_options','struct'));
  
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
