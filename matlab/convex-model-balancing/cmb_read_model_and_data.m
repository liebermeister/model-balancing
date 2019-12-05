function [network, q_info, data, c_init, kinetic_data, state_data] = cmb_read_model_and_data(cmb_options, file, options)

% [network, q_info, data, c_init, kinetic_data, state_data] = cmb_read_model_and_data(cmb_options, file, columns_mean, columns_std)
%
% file: either a single sbtab filename or a cell array of such filenames (files to be combined)
% options:      struct with (optional) fields
%      options.metabolite_table_id (default 'MetaboliteConcentration')
%      options.flux_table_id       (default 'Flux')
%      options.enzyme_table_id     (default 'EnzymeConcentration')
%      options.columns_mean:       cell array of column names for mean data values
%                                  (same columns names in metabolite, flux, and enzyme table!)
%      options.columns_std:        same, for std dev columns

eval(default('options','struct'));
  
% --------------------------------------------------------------
% Load data files -> network, kinetic_data, state_data
% --------------------------------------------------------------

result = sbtab_load_network_model_and_data(file, options);

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
