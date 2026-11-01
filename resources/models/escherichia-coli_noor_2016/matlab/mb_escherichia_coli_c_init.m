function [c_init, ind_e_off, v_init, Keq_init] = mb_escherichia_coli_c_init(model,mb_options)

% Set model name and input filenames

model_file = [mb_resourcedir '/models/escherichia-coli_noor_2016/escherichia-coli_noor_2016.tsv'];
network    = network_import_model(model_file, struct('load_quantity_table',0));

% set simple default concentrations

load([mb_basedir '/resources/models/escherichia-coli_noor_2016/escherichia-coli_example_reference_state']); % load 'c', 'u', 'v', 'Keq'

v_init = v;
c_init = c.emc4cm;
ind_e_off = network_find_reactions_by_name(network,{'FBP_R00762'});
Keq_init = Keq;

% c_init               = 0.1 * ones(length(network.metabolites),1);
%ind_high             = label_names({'D_Glucose'},network.metabolites);
%ind_high             = label_names({'D_Glucose','Ubiquinol','Acetyl_CoA'},network.metabolites);
% ind_low              = label_names({'ADP','Orthophosphate','NADH'},network.metabolites);
% ind_very_low         = label_names({'CO2','Ubiquinone'},network.metabolites);
%c_init(ind_high)     = 20;
% c_init(ind_low)      = 0.01;
% c_init(ind_very_low) = 0.001;

% ---------------------------------------------------------------------------------------
% Determine reference concentrations by enzyme cost minimization
% ---------------------------------------------------------------------------------------

if 0,
  
  data_dir                 = [ecm_resourcedir filesep 'models' filesep 'escherichia-coli_noor_2016'];
  filename_model           = [data_dir filesep 'escherichia-coli_noor_2016_ECM_Model.tsv'];
  filename_validation_data = [data_dir filesep 'escherichia-coli_noor_2016_ECM_ValidationData.tsv'];
  result_dir               = tempdir; 
  
  [network,v,c_data,u_data, conc_min, conc_max, met_fix, conc_fix,positions, enzyme_cost_weights, warnings] = load_model_and_data_sbtab(filename_model, filename_validation_data);
  
  ecm_options            = ecm_default_options(network, 'E. coli central carbon metabolism');
  ecm_options.c_data     = c_data;
  ecm_options.u_data     = u_data;
  ecm_options.ecm_scores = {'emc4cm'};
  ecm_options            = ecm_update_options(network, ecm_options);
  
  [c, u, u_cost, up, A_forward] = ecm_enzyme_cost_minimization(network, network.kinetics, v, ecm_options);
  Keq = network.kinetics.Keq;
  
  save([mb_basedir '/resources/models/escherichia-coli_noor_2016/escherichia-coli_example_reference_state'], 'c', 'u', 'v','Keq');

end