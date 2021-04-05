% ======================
% Model balancing
% ======================
% 
% Matlab files for model balancing
%
% o Main function: model_balancing.m
% o For usage examples, see m-files 'demo_cmb_artificial_data' and 'demo_cmb_experimental_data' 
% o For algorithm options, see 'cmb_default_options'
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
% kinetic_model.tsv    SBtab - Metabolic network and kinetic constants (optionally: include kinetic rate law formulae)
% metabolic_states.tsv SBtab - Predicted metabolite concentrations, enzyme concentrations, Gibbs free energies, fluxes (for all states)
% report.txt           text  - Number of data points / of fitted variables (by category), and calculation time
% options.tsv          SBtab - Options
% results.mat          .mat  - Results
