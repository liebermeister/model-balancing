% Tests for model balancing with a small model and artificial data

setup_test_path;
rng(1);
models  = [mb_resourcedir '/models/branch_point_model/'];
c_init  = [1, 1, 0.1, 0.5]';

mb_options                        = mb_default_options(1);
mb_options.ns                     = 2;
mb_options.prior_variant          = 'broad_prior';
mb_options.use_kinetic_data       = 'all';
mb_options.use_artificial_noise   = 0;
mb_options.use_kinetic_data_noise = 0;
mb_options.enzyme_score_alpha     = 0.5;
mb_options.initial_values_variant = 'average_sample';
mb_options.show_graphics          = 0;
mb_options.save_results           = 0;
mb_options.save_graphics          = 0;
mb_options.verbose                = 0;
mb_options.display                = 0;

[network, bounds, prior, q_info, data, true_values] = mb_model_artificial_data( ...
    [models 'branch_point_model.tsv'], mb_options, c_init, ...
    [models 'branch_point_model_Position.tsv'], [models 'branch_point_model_ConcentrationConstraint.tsv']);

filenames = mb_filenames('branch_point_model', 'test', 'test', tempdir, [models 'branch_point_model.tsv'], 1);

%% kinetic constants: conversion to the vector q and back
q = mb_kinetics_to_q(network, mb_options, q_info);
kinetics = mb_q_to_kinetics(q, network, mb_options, q_info);
assert(max(abs(kinetics.KV - network.kinetics.KV)) < 1e-8);
assert(max(abs(kinetics.Keq - network.kinetics.Keq)) < 1e-8);

%% model balancing gives finite results with the expected dimensions
optimal = model_balancing(filenames, mb_options, network, q_info, prior, bounds, data, true_values);
[nr, nm] = network_numbers(network);
assert(isequal(size(optimal.C), [nm, 2]));
assert(all(isfinite(optimal.C(:))));
