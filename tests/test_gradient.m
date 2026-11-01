% Test: the analytic gradient agrees with central finite differences
% (small models with regulation, a zero flux, a closed cycle, and different options;
% models from the folder matlab-python-comparison/test_model)

setup_test_path;
tm = fullfile(fileparts(mfilename('fullpath')), '..', 'matlab-python-comparison', 'test_model');
fn = @(name) fullfile(tm, ['test_model_' name '.tsv']);

test_models(1).files   = {fn('network'), fn('kinetic_data'), fn('state_data')};
test_models(1).columns = {'S1_Mean','S2_Mean','S3_Mean'};
test_models(2).files   = {fn('cycle_network'), fn('cycle_kinetic_data'), fn('cycle_state_data')};
test_models(2).columns = {'S1_Mean','S2_Mean'};

% options: alpha, beta, use_crosscovariances
variants = [0.5 0 0; 1 0 0; 0 0 0; 0.5 1 0; 0.5 0 1];

for im = 1:length(test_models)
  for iv = 1:size(variants,1)
    rng(1);
    o = mb_default_options;
    o.run = 'test';
    o.enzyme_score_alpha   = variants(iv,1);
    o.beta_ln_c_over_km    = variants(iv,2);
    o.use_crosscovariances = variants(iv,3);
    o.save_results = 0; o.save_graphics = 0; o.show_graphics = 0;
    file_options = struct('match_data_by','ModelElementId','columns_mean',{test_models(im).columns},'columns_std',{{}});
    [network,q_info,data,c_init] = mb_load_and_convert_model_and_data(o, test_models(im).files, file_options);
    [bounds,prior] = mb_model_and_data('test',network,data,q_info,c_init,o);
    pre = mb_prepare_posterior(prior,data,o,q_info);
    pp  = mb_make_pp(network);
    V   = data.V.mean;

    % a feasible point that is not an optimum: the polytope center
    o_start = o; o_start.enzyme_score_alpha = 0.5;
    [~,~,init] = mb_estimation(network,q_info,bounds,data,prior,pre,pp,V,o_start);
    [nm,ns] = size(init.X);
    y = mb_qX_to_y(init.q,init.X,nm,ns);

    [~,gradient] = mb_log_posterior(y,pp,pre,V,o,q_info);
    numerical = zeros(size(y));
    delta = 1e-5;
    for it = 1:length(y)
      dy = zeros(size(y)); dy(it) = delta;
      numerical(it) = [mb_log_posterior(y+dy,pp,pre,V,o,q_info) - mb_log_posterior(y-dy,pp,pre,V,o,q_info)] / (2*delta);
    end
    relative_error = max(abs(gradient-numerical)) / (1 + max(abs(numerical)));
    fprintf('model %d, variant %d: largest relative gradient error %.2e\n', im, iv, relative_error);
    assert(relative_error < 1e-5);
  end
end
