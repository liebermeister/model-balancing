function [bounds, prior, init] = mb_model_and_data(model, network, data, q_info, c_init, mb_options, conc_min, conc_max)

% [bounds, prior, init] = mb_model_and_data(model, network, data, q_info, c_init, mb_options, conc_min, conc_max)
%
% Prepare data structures describing a Convex Model Balancing problem 
%
% Input variables
%   model        (string) model name  
%   network      struct describing metabolic network (format as in Metabolic Network Toolbox)
%   data         struct describing data used in the optimality problem
%   q_info       struct describing the dependencies between model variables
%   c_init       matrix; initial guess of (non-logarithmic) metabolite concentrations
%   mb_options  struct containing options 
%   conc_min, conc_max vectors of minimal and maximal metabolite concentrations
%
% Output variables
%   bounds    struct describing bounds in the optimality problem (details see mb_make_bounds)
%   prior     struct describing priors in the optimality problem (details see mb_make_prior)
%   init      struct containing the initial model variables (fields init.q and init.X)
% 
% For problems with artificial data, use mb_model_artificial_data instead

eval(default('conc_min','[]','conc_max','[]'));
  
% Number of metabolic states
  
ns = size(data.X.mean,2);
mb_options.ns = ns;

% Bounds

bounds = mb_make_bounds(network,q_info,mb_options, conc_min, conc_max);

% Prior

prior = mb_make_prior(network,q_info,mb_options);

% Initial state for optimisation

init.q = mb_kinetics_to_q(network, mb_options, q_info);
init.X = log(c_init);
