function [network, data, kinetic_data] = cmb_read_model_and_data(data_set,version)

% Wrapper function for older function 'es_load_data';
  
par.model_version = version;
  
[~,~,~, condition, data, v_sign, n_exp, kinetic_data, parameter_prior, network, network_CoSplit, network_CoHid, ind_biomass, ind_biomass_production,N,W,nm,nr,ind_ext,kin_data] = es_load_data(data_set,par);
