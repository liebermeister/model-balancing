function [network, data, kinetic_data] = cmb_read_prepared_model_and_data(data_set,version)

% Wrapper function for older function 'es_load_data' - currently not working!
%
% Data sets supported:
%  'ycm_diauxic',
%  'ycm_oxygen',
%  'eccm_ishii_growth',
%  'eccm_old_ishii_growth',
%  'eccm_ishii_all',
%  'bscm_bigexp_MG10',
%  'bscm_bigexp',
%  'bscm_small_bigexp',

eval(default('version','[]'));
options.model_version = version;
  
[~,~,~, condition, data, v_sign, n_exp, kinetic_data, parameter_prior, network, network_CoSplit, network_CoHid, ind_biomass, ind_biomass_production,N,W,nm,nr,ind_ext,kin_data] = es_load_data(data_set,options);
