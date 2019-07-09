function state_data_files = e_coli_glucose_data_filenames(cmb_options)

% State data files for one mebtabolic state:

DATA_DIR  = [cmb_basedir '/resources/escherichia_coli/'];

state_data_files.conditions = {'Glucose'};

state_data_files.metabolite.file         = [DATA_DIR 'data-omics/growth_on_glucose/concentration_eco_glucose.tsv'];
state_data_files.metabolite.type         = 'metabolite_concentration';
state_data_files.metabolite.columns_mean = {'Value'};
state_data_files.metabolite.columns_std  = {};

state_data_files.enzyme.file             = [DATA_DIR 'data-omics/growth_on_glucose/enzyme_data.tsv'];
state_data_files.enzyme.type             = 'enzyme_concentration';
state_data_files.enzyme.columns_mean     = {'data'};
state_data_files.enzyme.columns_std      = {};

% haverkorn data, like in ECM paper
state_data_files.flux.file               = [DATA_DIR 'data-omics/growth_on_glucose/fluxes_eco_Chen_2010_aerobic.tsv'];
state_data_files.flux.type               = 'reaction_flux';
state_data_files.flux.columns_mean       = {'Value'};
state_data_files.flux.columns_std        = {'Std'};

% Alternative (leads to problems in calculations .. fix!)
% state_data_files.flux.file               = [DATA_DIR 'data-omics/growth_on_glucose/fluxes_eco_Haverkorn_2011_glucose_KEGG.tsv'];
% state_data_files.flux.columns_std        = {'Quantile95'};
