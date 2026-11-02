% Tests for the package setup: version number, default options, prior table

setup_test_path;

%% version number
assert(~isempty(regexp(model_balancing_version, '^\d+\.\d+\.\d+$', 'once')));

%% the prior table of this package exists and can be read
assert(exist(mb_prior_file, 'file') == 2);
parameter_prior = parameter_balancing_prior([], mb_prior_file, 0);
assert(any(strcmp(parameter_prior.Symbol, 'KM')));

%% default options
mb_options = mb_default_options;
assert(isstruct(mb_options));
