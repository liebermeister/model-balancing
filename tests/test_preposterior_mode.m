% Test: the initial value variant 'preposterior_mode' reaches the same optimum as 'polytope center'
% (E. coli model, kinetic data with alpha = 0.5 and alpha = 1; takes a few minutes)
% The constrained preposterior mode lies on the boundary of the feasible region. Without
% moving it into the interior, fmincon stops at the initial point.

setup_test_path;
m     = [mb_resourcedir '/models/e_coli_noor_2016/'];
files = {[m 'e_coli_noor_2016.tsv'], [m 'e_coli_ccm_kinetic_data/e_coli_noor_2016_kinetic_data_modified.tsv'], ...
         [m 'e_coli_ccm_state_data_COPY_FOR_INFO.tsv']};
file_options.match_data_by = 'ModelElementId';
file_options.columns_mean  = {'S_1_Mean'};
file_options.columns_std   = {};
file_options.constraint_sbtab_file = [m 'e_coli_noor_2016_ConcentrationConstraint.tsv'];

for alpha = [0.5, 1]
  value = [nan, nan];
  variants = {'polytope center', 'preposterior_mode'};
  for it = 1:2
    rng(1);
    o = mb_default_options;
    o.run = 'test';
    o.enzyme_score_alpha = alpha;
    o.initial_values_variant = variants{it};
    o.save_results = 0; o.save_graphics = 0; o.show_graphics = 0;
    [network,q_info,data,c_init,~,~,conc_min,conc_max] = mb_load_and_convert_model_and_data(o,files,file_options);
    [bounds,prior] = mb_model_and_data('test',network,data,q_info,c_init,o,conc_min,conc_max);
    pre = mb_prepare_posterior(prior,data,o,q_info);
    pp  = mb_make_pp(network);
    opt = mb_estimation(network,q_info,bounds,data,prior,pre,pp,data.V.mean,o);
    [nr,nm] = network_numbers(network);
    value(it) = mb_log_posterior(mb_qX_to_y(opt.q,opt.X,nm,size(data.V.mean,2)),pp,pre,data.V.mean,o,q_info);
  end
  fprintf('alpha %g: polytope center %.4f, preposterior_mode %.4f\n',alpha,value(1),value(2));
  assert(abs(value(1) - value(2)) < 0.01);
end
