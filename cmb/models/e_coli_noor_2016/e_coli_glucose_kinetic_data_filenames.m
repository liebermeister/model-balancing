function [model_file, kinetic_data_file] = e_coli_glucose_kinetic_data_filenames(cmb_options)

% Kinetic data files

DATA_DIR  = [cmb_basedir '/resources/escherichia_coli/'];

model_file = [ DATA_DIR '/network/ecoli_noor_2016.tsv'];

switch cmb_options.source_of_kinetic_data, 
  
  case 'original',    
    %% original data in one file, written by script e_coli_balance_parameters.m

    display('Using in-vitro kinetic constant data');    
    kinetic_data_file = [DATA_DIR 'data-kinetic/e_coli_balanced_parameters_preprocessed_data.tsv'];

    %% ALTERNATIVE: from original data files
    %% kinetic_data_file = {[DATA_DIR 'data-kinetic/ecoli_ccm_equilibrium_constants_from_ECM_model.tsv'], ...
    %%                      [DATA_DIR 'data-kinetic/kinetic_constants_brenda_eco.tsv'], ...
    %%                      [DATA_DIR 'data-kinetic/kinetic_constants_flamholz_several_organisms.tsv'], ...
    %%                      [DATA_DIR 'data-kinetic/eco_gapdh_data.tsv'], ...
    %%                      [DATA_DIR 'data-kinetic/eco_pts_data.tsv']};

  case 'balanced',
    %%  balanced by script e_coli_balance_parameters.m

    display('Using balanced kinetic constants from E. coli model as surrogate kinetic data');
    kinetic_data_file = [DATA_DIR 'data-kinetic/e_coli_balanced_parameters.tsv'];

    %% OUTDATED: use balanced parameters from Noor et al 2016
    %%
    %% dum = sbtab_document_load([DATA_DIR 'data-kinetic/ecoli_ccm_ProteinComposition_Haverkorn_ModelState.tsv']);   
    %% kinetic_data_file = sbtab_document_get_table(dum,'Parameter');

  otherwise,
    error('');
end
