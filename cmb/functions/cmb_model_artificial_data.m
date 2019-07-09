function [filenames, network, q_info, prior, bounds, data, true] = cmb_model_artificial_data(model_name, network_sbml, cmb_options, c_init)

% set filenames

DIR                      = [cmb_basedir '/models/' model_name '/'];
  
filenames.model_name     = model_name;
filenames.network_sbml   = network_sbml;
filenames.run_dir        = [ DIR  'simulations/' cmb_options.run ];
filenames.graphics_dir   = [ DIR  'simulations/' cmb_options.run '/ps-files/' ];
filenames.data           = [ DIR  'simulations/' cmb_options.run '/data/'     ];
filenames.states_out     = [ DIR  'simulations/' cmb_options.run '/data/'];
filenames.parameters_out = [ DIR  'simulations/' cmb_options.run '/data/parameters.tsv'];
filenames.model_out      = [ DIR  'simulations/' cmb_options.run '/data/model.tsv'];
filenames.report         = [ DIR  'simulations/' cmb_options.run '/data/report.txt'];
filenames.options_file   = [ DIR  'simulations/' cmb_options.run '/data/options'];
filenames.result_file    = [ DIR  'simulations/' cmb_options.run '/data/results'];

[~,~] = mkdir(filenames.graphics_dir);
[~,~] = mkdir(filenames.data);

% load network 

network = network_sbml_import(filenames.network_sbml);

% parameter structure (depends on cmb_options.parameterisation)

q_info = cmb_define_parameterisation(network, cmb_options); 

% generate model kinetics and artificial data

[kinetics, prior, bounds, data, true] = cmb_generate_artificial_data(network, cmb_options, q_info, c_init);

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
