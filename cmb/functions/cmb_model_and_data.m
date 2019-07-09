function [ns, bounds, prior, init, filenames] = cmb_model_and_data(model, network, data, q_info, c_init, cmb_options, cmb_DIR)

% Number of metabolic states
  
ns = size(data.X.mean,2);
cmb_options.ns = ns;

% Bounds

bounds = cmb_make_bounds(network,q_info,cmb_options);

% Prior

prior = cmb_make_prior(network,q_info,cmb_options);

% Initial state for optimisation

init.q = cmb_kinetics_to_q(network, cmb_options, q_info);
init.X = log(c_init);

% Filenames for convex estimation runs

filenames.model_name     = model;
filenames.graphics_dir   = [ cmb_DIR '/simulations/' cmb_options.run '/ps-files/' ];
filenames.data           = [ cmb_DIR '/simulations/' cmb_options.run '/data/' ];
filenames.states_out     = [ cmb_DIR '/simulations/' cmb_options.run '/data/'];
filenames.parameters_out = [ cmb_DIR '/simulations/' cmb_options.run '/data/parameters.tsv'];
filenames.model_out      = [ cmb_DIR '/simulations/' cmb_options.run '/data/model.tsv'];
filenames.report         = [ cmb_DIR '/simulations/' cmb_options.run '/data/report.txt'];
filenames.options_file   = [ cmb_DIR '/simulations/' cmb_options.run '/data/options'];
filenames.result_file    = [ cmb_DIR '/simulations/' cmb_options.run '/data/results'];

[~,~] = mkdir(filenames.graphics_dir);
[~,~] = mkdir(filenames.data);
