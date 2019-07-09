function state_data_files = e_coli_data_filenames(cmb_options)

% State data files for one metabolic state:

DATA_DIR  = [cmb_basedir '/resources/escherichia_coli'];

% THERE ARE MORE glucose conditions - change this!!
state_data_files.samples = {'Glucose','Glycerol','Acetate'}';

state_data_files.metabolite.file         = [DATA_DIR '/data-omics/davidi_2016/Davidi_2016_metabolite_concentrations_uM.tsv'];
state_data_files.metabolite.type         = 'metabolite_concentration';
state_data_files.metabolite.columns_mean = {'>Condition:Glucose', '>Condition:Glycerol', '>Condition:Acetate'};
state_data_files.metabolite.columns_std  = {};

state_data_files.enzyme.file             = [DATA_DIR '/data-omics/davidi_2016/Davidi_2016_protein_abundance_mmolgCDW_KEGG_IDs.tsv'];
state_data_files.enzyme.type             = 'enzyme_concentration';
state_data_files.enzyme.columns_mean     = {'GLC_CHEM_mu_0_11_V','GLYC_BATCH_mu_0_47_S', 'ACE_BATCH_mu_0_3_S'};
state_data_files.enzyme.columns_std      = {};
                              
state_data_files.flux.file               = [DATA_DIR '/data-omics/davidi_2016/Davidi_2016_flux_mmolgCDW_s_KEGG_IDs.tsv'];
state_data_files.flux.type               = 'reaction_flux';
state_data_files.flux.columns_mean       = {'GLC_CHEM_mu_0_11_V','GLYC_BATCH_mu_0_47_S', 'ACE_BATCH_mu_0_3_S'};
state_data_files.flux.columns_std        = {};
