% Tests for the handling of missing data, kinetic data errors, and networks with closed cycles

setup_test_path;

cmb_options = cmb_default_options(1);
cmb_options.ns = 2;

% network with a closed cycle (A -> B -> C and A -> C): one Wegscheider condition
network = network_build_from_sum_formulae_list({'Xa <=> A','A <=> B','B <=> C','A <=> C','C <=> Xb'});
network.kinetics = set_kinetics(network,'cs');
[nr, nm] = network_numbers(network);
q_info = cmb_define_parameterisation(network, cmb_options);
bounds = cmb_make_bounds(network, q_info, cmb_options);
prior  = cmb_make_prior(network, q_info, cmb_options);

% data: no kinetic data, all metabolites and enzymes measured except for one of each
data.X.mean   = log(ones(nm,2));
data.X.std    = 0.2 * ones(nm,2);
data.lnE.mean = log(ones(nr,2));
data.lnE.std  = 0.2 * ones(nr,2);
data.V.mean   = ones(nr,2);
data.V.std    = nan(nr,2);
data.qall.mean = nan(q_info.qall.number,1);
data.qall.std  = nan(q_info.qall.number,1);
data.X.mean(1,2)   = nan;
data.lnE.mean(2,1) = nan;

preposterior = cmb_prepare_posterior(prior, data, cmb_options, q_info);

%% missing metabolite data: mean and std of the prior
assert(abs(preposterior.X.mean(1,2) - prior.X.mean(1,2)) < 1e-12);
assert(abs(preposterior.X.std(1,2)  - prior.X.std(1,2)) < 1e-12);

%% missing enzyme data: mean and std of the prior
assert(abs(preposterior.lnE.mean(2,1) - prior.lnE.mean(2,1)) < 1e-12);
assert(abs(preposterior.lnE.std(2,1)  - prior.lnE.std(2,1)) < 1e-12);

%% elements with data are narrower than the prior
assert(preposterior.X.std(1,1) < 0.2);
assert(preposterior.lnE.std(2,2) < 0.2);

%% closed cycle: fewer independent than all equilibrium constants
assert(length(q_info.q.index.Keq) == nr - 1);

%% closed cycle: bounds have the length of the vector qall
assert(length(bounds.q_all_min) == q_info.qall.number);
assert(length(bounds.q_all_max) == q_info.qall.number);

%% closed cycle: all Keq values obey the Wegscheider condition
q = randn(q_info.q.number,1);
kinetics = cmb_q_to_kinetics(q, network, cmb_options, q_info);
null_space = null(full(network.N));
assert(max(abs(null_space' * log(kinetics.Keq))) < 1e-10);

%% options: the upper bound on KA comes from the row of KA in the prior table
parameter_prior = parameter_balancing_prior([], cmb_prior_file, 0);
upper_bound = cell_string2num(parameter_prior.UpperBound);
assert(cmb_options.quantities.KA.max == upper_bound(strcmp(parameter_prior.Symbol,'KA')));

%% kinetic data: the standard deviation of log values is positive and given by the geometric std
kinetic_data = kinetic_data_load([], parameter_balancing_prior([], cmb_prior_file, 0), network, [], struct('enforce_ranges',0));
kinetic_data.Keq.median(1) = 2;
kinetic_data.Keq.std_ln(1) = log(1.5);
state_data.metabolite_data = struct('Mean',ones(nm,2),'Std',[]);
state_data.enzyme_data     = struct('Mean',ones(nr,2),'Std',[]);
state_data.flux_data       = struct('Mean',ones(nr,2),'Std',[]);
d = state_data_to_data(kinetic_data, state_data, cmb_options, network);
assert(abs(d.qall.std(q_info.qall.index.Keq(1)) - log(1.5)) < 1e-12);
