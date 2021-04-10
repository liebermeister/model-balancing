function c_init = e_coli_c_init(model,cmb_options)

% Set model name and input filenames

model_file = [cmb_resourcedir '/models/e_coli_noor_2016/e_coli_noor_2016.tsv'];
network    = network_import_model(model_file, struct('load_quantity_table',0));

% set simple default concentrations

c_init       = 0.01 * ones(length(network.metabolites),1);
ind_glucose  = label_names('D_Glucose',network.metabolites);
ind_products = label_names({'CO2','Ubiquinone'},network.metabolites);
c_init(ind_glucose)  = 10;
c_init(ind_products) = 0.00001;
