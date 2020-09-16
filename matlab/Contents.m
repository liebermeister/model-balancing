% ======================
% Convex model balancing
% ======================
% 
% Matlab files for convex model balancing
% Main function: convex_model_balancing.m
%
% For an example, see 'demo_cmb_artificial_data.m'
% 
% ------------
% Input data
% ------------
% 
% Model
% 
% Model structures must be given as SBtab or SBML files
% 
% Kinetic data
% 
% Kinetic data must be given as SBtab files
% o Metabolites must be annotated with KEGG Compound IDs
% o Reactions and Enzymes must be annotated with KEGG Reaction IDs
% 
% Metabolic state data
% 
% Metabolite, flux, and enzyme data must be given as SBtab files
% o Metabolites must be annotated with KEGG Compound IDs
% o Reactions and Enzymes must be annotated with KEGG Reaction IDs
% 
% The columns containing the data can carry arbitrary names 
% (which must be declared within matlab, for parsing the files)
% 
% ------------
% Result files
% ------------
% 
% kinetic_model.tsv    SBtab - Metabolic network and kinetic constants (optionally: include kinetic rate law formulae)
% metabolic_states.tsv SBtab - Predicted metabolite concentrations, enzyme concentrations, Gibbs free energies, fluxes (for all states)
% report.txt           Number of data points / of fitted variables (by category), and calculation time
% 
% options.mat      .mat  Options
% results.mat      .mat  Results
