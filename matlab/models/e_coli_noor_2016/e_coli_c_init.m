function c_init = e_coli_c_init(model,cmb_options)

% --------------------------------------------------------------
% Set model name and input filenames

model_file = [cmb_resourcedir '/data/data-organisms/escherichia_coli/network/ecoli_noor_2016.tsv'];

network = network_import_model(model_file);

% set simple default concentrations

ind_glucose = label_names('D_Glucose',network.metabolites);
c_init                         = 0.1 * ones(length(network.metabolites),1);
c_init(find(network.external)) = 0.001;
c_init(ind_glucose)            = 10;
