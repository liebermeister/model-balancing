function [network, bounds, prior, q_info, data, true, kinetic_data, state_data, conc_min, conc_max] = cmb_model_artificial_data(network_file, cmb_options, c_init, position_sbtab_file, constraint_sbtab_file)

% [network, bounds, prior, q_info, data, true, kinetic_data, state_data] = cmb_model_artificial_data(network_file, cmb_options, c_init, position_sbtab_file, constraint_sbtab_file)
%
% Generate input data for a Convex Model Balancing problem with artifical data
% 
% Input variables
%   network_file filename of SBML or SBtab network model
%                an SBtab file must contain tables "Reaction", "Compound" (and optionally "Position")
%   cmb_options  struct containing options
%   .use_artificial_noise    = 0;                   % flag: generate artificial state data with noise
%   .use_kinetic_data_noise  = 1;                   % flag: generate artificial kinetic data with noise
%   .prior_variant           = 'original_prior';                original prior
%                              'broad_prior'                    broad prior around original prior means
%                              'broad_prior_around_zero'        broad prior around 0 values
%                              'prior_around_true_values'       prior around true values (only for artificial data)
%                              'broad_prior_around_true_values' broad prior around true values (only for artificial data)

%   c_init       (optional) matrix; initial guess of (non-logarithmic) metabolite concentrations; 
%                           also used for generating the artificial data 
%   position_sbtab_file (optional) filename of SBtab file (table "Position") for network layout 
%                                  (overrides an existing table "Position" in network_file)
%   constraint_sbtab_file (optional) filename of SBtab file (tables "Position") 
% 
% Output variables
%   network      struct describing metabolic network (format as in Metabolic Network Toolbox)
%   bounds       struct describing bounds in the optimality problem (details see cmb_make_bounds)
%   prior        struct describing priors in the optimality problem (details see cmb_make_prior)
%   q_info       struct describing the dependencies between model variables
%   data         struct describing data used in the optimality problem
%   true         struct describing true model variables (optional; only for models with artificial data)
%   kinetic_data
%   state_data
%
% For CMB problems with "real" data, use cmb_model_and_data instead
  
eval(default('c_init','[]','position_sbtab_file','[]'));

% load network 

if cmb_options.verbose,
  display(sprintf('Reading file %s', network_file));
end

network = network_import(network_file);

if length(position_sbtab_file),
  network = netgraph_read_positions(network, position_sbtab_file);
end

% concentration bounds
conc_min = [];
conc_max = [];

if length(constraint_sbtab_file),
  constraint_sbtab = sbtab_document_load_from_one(constraint_sbtab_file);
  conc_min  = sbtab_table_get_column(constraint_sbtab.tables.ConcentrationConstraint,'Min',1);
  conc_max  = sbtab_table_get_column(constraint_sbtab.tables.ConcentrationConstraint,'Max',1);
end

% parameter structure (depends on cmb_options.parameterisation)

q_info = cmb_define_parameterisation(network, cmb_options);

% generate model kinetics and artificial data

[kinetics, prior, bounds, data, true, kinetic_data, state_data] = cmb_generate_artificial_data(network, cmb_options, q_info, c_init, conc_min, conc_max);

network.kinetics = kinetics;

% possibly, use different kinetic constants prior (for reconstruction) than prior used to generate the model kinetics

switch cmb_options.prior_variant,
  
  case 'original_prior',
  
  case 'prior_around_true_values',
    prior.q.mean = true.q;
    prior.qdep_pseudo.mean = q_info.M_q_to_qdep * true.q; 
    prior.X.mean = true.X;
    prior.lnE.mean = log(true.E);

  case 'broad_prior',
    prior.q.std   = log(10^10)*ones(size(prior.q.std));
    prior.qdep_pseudo.std  = log(10^10)*ones(size(prior.qdep_pseudo.std));
    prior.X.std   = log(10^10)*ones(size(prior.X.std));
    prior.lnE.std = log(10^10)*ones(size(prior.lnE.std));
  
  case 'broad_prior_around_zero',
    prior.q.mean = zeros(size(prior.q.std));
    prior.q.std  = log(10^10)*ones(size(prior.q.std));
    prior.qdep_pseudo.mean = zeros(size(prior.qdep_pseudo.std)); 
    prior.qdep_pseudo.std  = log(10^10)*ones(size(prior.qdep_pseudo.std));
    prior.X.mean = 0 * prior.X.mean;
    prior.X.std  = log(10^10)*ones(size(prior.X.std));
    prior.lnE.mean = 0 * prior.lnE.mean;
    prior.lnE.std  = log(10^10)*ones(size(prior.lnE.std));
  
  case 'broad_prior_around_true_values',
    prior.q.mean = true.q;
    prior.q.std  = log(10^10)*ones(size(prior.q.std));
    prior.qdep_pseudo.mean = q_info.M_q_to_qdep * true.q; 
    prior.qdep_pseudo.std  = log(10^10)*ones(size(prior.qdep_pseudo.std));
    prior.X.mean = true.X;
    prior.X.std  = log(10^10)*ones(size(prior.X.std));
    prior.lnE.mean = true.lnE;
    prior.lnE.std  = log(10^10)*ones(size(prior.lnE.std));

  otherwise, error('');
end

prior.q.prec           = diag(1./prior.q.std.^2);
prior.qdep_pseudo.prec = diag(1./prior.qdep_pseudo.std.^2);
