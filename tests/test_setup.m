% Tests for the package setup: version number, default options, prior table

setup_test_path;

%% version number
assert(~isempty(regexp(model_balancing_version, '^\d+\.\d+\.\d+$', 'once')));

%% the prior table of this package exists and can be read
assert(exist(cmb_prior_file, 'file') == 2);
parameter_prior = parameter_balancing_prior([], cmb_prior_file, 0);
assert(any(strcmp(parameter_prior.Symbol, 'KM')));

%% default options
cmb_options = cmb_default_options;
assert(isstruct(cmb_options));
assert(strcmp(cmb_options.parameterisation, 'Keq_KV_KM_KA_KI'));
