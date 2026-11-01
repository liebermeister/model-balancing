% ======================
% Model balancing
% ======================
% 
% Matlab files for model balancing
%
% o Main function: model_balancing.m
% o For usage examples, see m-files 'demo_mb_artificial_data' and 'demo_mb_experimental_data' 
% o For algorithm options, see 'mb_default_options'
% 
% ------------
% Input data
% ------------
% 
% Model:
%
%  Model structures are read from SBtab or SBML files
% 
% Kinetic data:
% 
%  Kinetic data are read from SBtab files
%  o Metabolites must carry the same IDs as in the model (or KEGG Compound IDs)
%  o Reactions and Enzymes must carry the same IDs as in the model (or KEGG Reaction IDs)
% 
% Metabolic state data:
% 
%  Metabolite, flux, and enzyme data are read from SBtab files
%  o Metabolites must carry the same IDs as in the model (or KEGG Compound IDs)
%  o Reactions and Enzymes must carry the same IDs as in the model (or KEGG Reaction IDs)
%  
%  The columns containing the data can carry arbitrary names 
%    (which must be declared within matlab, for parsing the files)
% 
% ------------
% Output files
% ------------
% 
% Output files are written to [result_dir]/[model]/[run]/data/ (see mb_filenames):
%
% balanced_model.tsv   SBtab - Metabolic network and balanced kinetic constants
% balanced_states.tsv  SBtab - Balanced metabolite concentrations, enzyme concentrations, Gibbs free energies, fluxes (for all states)
% data_kinetic.tsv     SBtab - Kinetic data used
% data_states.tsv      SBtab - Metabolic state data used
% true_model.tsv       SBtab - True model (only for artificial data)
% true_states.tsv      SBtab - True metabolic states (only for artificial data)
% report.txt           text  - Number of data points / of fitted variables (by category), and calculation time
% options.tsv          SBtab - Options
% results.mat          .mat  - Results
%
% Graphics files are written to [result_dir]/[model]/[run]/ps-files/
%
% ---------
% Functions
% ---------
%
% Main function (in the folder above): model_balancing
%
% Problem setup and options
%   mb_default_options                  - Default options for model balancing
%   mb_filenames                        - Filenames for input and output files
%   mb_model_and_data                   - Prepare a model balancing problem (bounds, prior, initial values)
%   mb_model_artificial_data            - Prepare a model balancing problem with artificial data
%   mb_generate_artificial_data         - Generate artificial kinetic data and state data for a network
%   mb_make_bounds                      - Bounds on model variables
%   mb_make_prior                       - Prior distributions for model variables
%   mb_make_pp                          - Struct with network information for enzyme cost functions
%   mb_project_fluxes                   - Project fluxes (e.g., to complete missing flux values)
%   mb_prior_file                       - Filename of the default prior table
%   mb_basedir                          - Main directory of the model balancing repository
%   mb_resourcedir                      - Directory with example models and data
%
% Estimation
%   mb_estimation                       - Run the optimisation (called by model_balancing)
%   mb_prepare_posterior                - Preposterior distributions (prior and data combined)
%   mb_objective                        - Objective function for fmincon
%   mb_log_posterior                    - Log posterior score
%   mb_variables_to_log_posterior       - Log posterior score for given kinetic constants and concentrations
%
% Results and graphics
%   mb_display_scores                   - Display log posterior scores (true, initial, optimised values)
%   mb_statistics                       - Report about a model balancing run
%   mb_graphics                         - Scatter plots: balanced vs. data values and vs. true values
%   mb_diagnostic_network_graphics      - Diagnostic network graphics
%
% Subfolder mb-variables: conversion between variables
%   mb_define_parameterisation          - Define the parameterisation (independent and dependent kinetic constants)
%   mb_kinetics_to_q, mb_q_to_kinetics - Convert between kinetics struct and vector q
%   mb_qall_to_kinetics                 - Convert vector qall into kinetics struct
%   mb_q_to_qall, mb_qall_to_q         - Convert between vectors q and qall
%   mb_qX_to_y, mb_y_to_qX             - Convert between (q, X) and optimisation vector y
%   mb_qX_to_variables                  - Compute all model variables from q, X, and fluxes
%
% Subfolder mb-file-input-output: file input and output
%   mb_load_model_and_data              - Load model and data from SBML / SBtab files
%   mb_load_and_convert_model_and_data  - Load model and data from SBtab files and convert them for model balancing
%   mb_save_results                     - Save model balancing results as SBtab files
%   mb_read_prepared_model_and_data     - Wrapper for older function es_load_data (currently not working)
%
% Subfolder mb-cvxpy: export for the python implementation
%   cvxpy_problem_data_structure         - Convert a model balancing problem into a struct for JSON export
%   cvxpy_combine_prior_and_likelihood   - Combine prior and likelihood (Gaussian distributions on log scale)
