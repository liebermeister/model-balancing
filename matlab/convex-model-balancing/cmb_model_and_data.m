function [ns, bounds, prior, init, filenames] = cmb_model_and_data(model, network, data, q_info, c_init, cmb_options, out_DIR)

% [ns, bounds, prior, init, filenames] = cmb_model_and_data(model, network, data, q_info, c_init, cmb_options, out_DIR)
%
% Prepare data structures describing a Convex Model Balancing problem 
%
% Input variables
%   model        (string) model name  
%   network      struct describing metabolic network (format as in Metabolic Network Toolbox)
%   data         struct describing data used in the optimality problem
%   q_info       struct describing the dependencies between model variables
%   c_init       matrix; initial guess of (non-logarithmic) metabolite concentrations
%   cmb_options  struct containing options 
%   out_DIR      Name of output directory (simple string, no directory path)    
%
% Output variables
%   ns        number of states (integer) 
%   bounds    struct describing bounds in the optimality problem (details see cmb_make_bounds)
%   prior     struct describing priors in the optimality problem (details see cmb_make_prior)
%   init      struct containing the initial model variables (fields init.q and init.X)
%   filenames struct describing filenames used (constructed using [out_DIR] and [cmb_options.run])
% 
% For problems with artificial data, use cmb_model_artificial_data instead

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
filenames.graphics_dir   = [ out_DIR '/simulations/' cmb_options.run '/ps-files/' ];
filenames.data           = [ out_DIR '/simulations/' cmb_options.run '/data/' ];
filenames.states_out     = [ out_DIR '/simulations/' cmb_options.run '/data/'];
filenames.parameters_out = [ out_DIR '/simulations/' cmb_options.run '/data/parameters.tsv'];
filenames.model_out      = [ out_DIR '/simulations/' cmb_options.run '/data/model.tsv'];
filenames.report         = [ out_DIR '/simulations/' cmb_options.run '/data/report.txt'];
filenames.options_file   = [ out_DIR '/simulations/' cmb_options.run '/data/options'];
filenames.result_file    = [ out_DIR '/simulations/' cmb_options.run '/data/results'];

[~,~] = mkdir(filenames.graphics_dir);
[~,~] = mkdir(filenames.data);
