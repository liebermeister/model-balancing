function [network, bounds, prior, q_info, data, true, kinetic_data, state_data] = cmb_model_artificial_data(network_file, cmb_options, c_init, position_sbtab_file)

% [network, bounds, prior, q_info, data, true, kinetic_data] = cmb_model_artificial_data(network_file, cmb_options, c_init, position_sbtab_file)
%
% Generate input data for a Convex Model Balancing problem with artifical data
% 
% Input variables
%   network_file filename of SBML or SBtab network model 
%   cmb_options  struct containing options
%   c_init       (optional) matrix; initial guess of (non-logarithmic) metabolite concentrations; 
%                           also used for generating the artificial data 
%   position_sbtab_file (optional) filename of SBtab file (table "Position") for network layout
% 
% Output variables
%   network      struct describing metabolic network (format as in Metabolic Network Toolbox)
%   q_info       struct describing the dependencies between model variables
%   prior        struct describing priors in the optimality problem (details see cmb_make_prior)
%   bounds       struct describing bounds in the optimality problem (details see cmb_make_bounds)
%   data         struct describing data used in the optimality problem
%   true         struct describing true model variables (optional; only for models with artificial data)
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

% parameter structure (depends on cmb_options.parameterisation)

q_info = cmb_define_parameterisation(network, cmb_options);

% generate model kinetics and artificial data

[kinetics, prior, bounds, data, true, kinetic_data, state_data] = cmb_generate_artificial_data(network, cmb_options, q_info, c_init);

network.kinetics = kinetics;

% possibly, use different kinetic constants prior (for reconstruction) than prior used to generate the model kinetics

switch cmb_options.prior_variant,
  
  case 'original_prior',
  
  case 'prior_around_true_values',
    prior.q.mean = true.q;
    prior.X.mean = true.X;
    prior.E.mean = true.E;

  case 'broad_prior',
    prior.q.std  = log(10^10)*ones(size(prior.q.std));
  
  case  'broad_prior_around_zero',
    prior.q.mean = zeros(size(prior.q.std));
    prior.q.std  = log(10^10)*ones(size(prior.q.std));
  
  case 'broad_prior_around_true_values',
    prior.q.mean = true.q;
    prior.q.std  = log(10^10)*ones(size(prior.q.std));

  otherwise, error('');
end

prior.q.cov_inv = diag(1./prior.q.std.^2);
