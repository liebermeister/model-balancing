% Comparison of the Matlab and Python versions of model balancing - Matlab part
%
% For each test case, this script
%   - writes the JSON problem file with the Matlab code
%   - runs model balancing
%   - saves the optimum, enzyme levels, and scores (file *_matlab.json)
% The Python scripts compare_json.py and compare_objective.py use these files.
% Output directory: environment variable MB_COMPARISON_OUT (default: folder 'results' next to this file)
%
% To run everything, see run_comparison.sh

here = fileparts(mfilename('fullpath'));
outdir = getenv('MB_COMPARISON_OUT');
if isempty(outdir), outdir = fullfile(here,'results'); end
[~,~] = mkdir(outdir);

% add the packages to the matlab path (this script should not depend on startup.m)
base = getenv('MB_MATLAB_PACKAGES');
if isempty(base), base = fileparts(fileparts(here)); end
for package = {'parameter-balancing/parameter-balancing','enzyme-cost-minimization/enzyme-cost-minimization', ...
               'model-balancing/model-balancing','sbtab/sbtab','structural-thermokinetic-modelling/stm', ...
               'matlab-utils/matlab-utils','mnt/mnt'}
  p = strsplit(genpath(fullfile(base,package{1})),pathsep);
  addpath(strjoin(p(~contains(p,[filesep '.'])),pathsep));
end

res = mb_resourcedir;
tm  = fullfile(here,'test_model');
fn  = @(name) fullfile(tm,['test_model_' name '.tsv']);
s3  = {'S1_Mean','S2_Mean','S3_Mean'};
s3s = {'S1_Std','S2_Std','S3_Std'};

test_model = struct('files',{{fn('network'),fn('kinetic_data'),fn('state_data')}},'constraints',fn('constraints'),'mean',{s3},'std',{s3s});

cases = {};

c = test_model; c.name = 'test_model_alpha_0.5';          c.alpha = 0.5; c.use_kinetic_data = 'all';      cases{1} = c;
c = test_model; c.name = 'test_model_alpha_1';            c.alpha = 1;   c.use_kinetic_data = 'all';      cases{2} = c;
c = test_model; c.name = 'test_model_alpha_0';            c.alpha = 0;   c.use_kinetic_data = 'all';      cases{3} = c;
c = test_model; c.name = 'test_model_no_std';             c.alpha = 0.5; c.use_kinetic_data = 'all';      c.std = {}; cases{4} = c;
c = test_model; c.name = 'test_model_only_Keq';           c.alpha = 0.5; c.use_kinetic_data = 'only_Keq'; cases{5} = c;
c = test_model; c.name = 'test_model_no_kinetic_data';    c.alpha = 0.5; c.use_kinetic_data = 'none';     cases{6} = c;

c = struct('files',{{fn('cycle_network'),fn('cycle_kinetic_data'),fn('cycle_state_data')}},'constraints',[],'mean',{{'S1_Mean','S2_Mean'}},'std',{{'S1_Std','S2_Std'}});
c.name = 'test_model_cycle'; c.alpha = 0.5; c.use_kinetic_data = 'all'; cases{7} = c;

m = [res '/models/branch_point_model/'];
c = struct('files',{{[m 'data/artificial_network_true.tsv'],[m 'data/artificial_kinetic_data.tsv'],[m 'data/artificial_state_data.tsv']}}, ...
           'constraints',[m 'branch_point_model_ConcentrationConstraint.tsv'], ...
           'mean',{{'S1_Mean','S2_Mean','S3_Mean','S4_Mean','S5_Mean','S6_Mean'}},'std',{{'S1_Std','S2_Std','S3_Std','S4_Std','S5_Std','S6_Std'}});
c.name = 'branch_point_model'; c.alpha = 0.5; c.use_kinetic_data = 'all'; cases{8} = c;

m = [res '/models/escherichia-coli_noor_2016/'];
c = struct('files',{{[m 'escherichia-coli_noor_2016.tsv'],[m 'escherichia-coli_ccm_kinetic_data/escherichia-coli_noor_2016_kinetic_data_modified.tsv'],[m 'escherichia-coli_ccm_state_data_COPY_FOR_INFO.tsv']}}, ...
           'constraints',[m 'escherichia-coli_noor_2016_ConcentrationConstraint.tsv'],'mean',{{'S_1_Mean'}},'std',{{}});
c.name = 'escherichia-coli_noor_2016'; c.alpha = 0.5; c.use_kinetic_data = 'all'; cases{9} = c;

manifest = {};

for it = 1:length(cases)
  c = cases{it};
  fprintf('\n--- %s\n',c.name);

  clear file_options
  file_options.match_data_by = 'ModelElementId';
  file_options.columns_mean  = c.mean;
  file_options.columns_std   = c.std;
  if ~isempty(c.constraints), file_options.constraint_sbtab_file = c.constraints; end

  rng(1); % the default initial point is computed from random LP solutions
  o = mb_default_options;
  o.run = c.name;
  o.enzyme_score_alpha = c.alpha;
  o.use_kinetic_data = c.use_kinetic_data;
  o.save_results = 0; o.save_graphics = 0; o.show_graphics = 0;

  [network,q_info,data,c_init,~,~,conc_min,conc_max] = mb_load_and_convert_model_and_data(o,c.files,file_options);
  [bounds,prior] = mb_model_and_data(c.name,network,data,q_info,c_init,o,conc_min,conc_max);

  mb_save_json_file(cvxpy_problem_data_structure(network,q_info,prior,data,bounds,[],o),fullfile(outdir,[c.name '_problem.json']));

  pre = mb_prepare_posterior(prior,data,o,q_info);
  pp  = mb_make_pp(network);
  V   = data.V.mean;
  opt = mb_estimation(network,q_info,bounds,data,prior,pre,pp,V,o);

  [nr,nm,~,KM_i,KA_i,KI_i] = network_numbers(network);
  ns = size(V,2);
  y  = mb_qX_to_y(opt.q,opt.X,nm,ns);
  o2 = o; o2.beta_ln_c_over_km = 2;
  r  = mb_q_to_kinetics(opt.q,network,o,q_info);

  result = struct('log_posterior',mb_log_posterior(y,pp,pre,V,o,q_info), ...
                  'log_posterior_beta_2',mb_log_posterior(y,pp,pre,V,o2,q_info), ...
                  'alpha',c.alpha,'X',{num2cell_matrix(opt.X)},'E',{num2cell_matrix(opt.E)},'V',{num2cell_matrix(V)}, ...
                  'lnEmean',{num2cell_matrix(pre.lnE.mean)},'lnEstd',{num2cell_matrix(pre.lnE.std)}, ...
                  'Xmean',{num2cell_matrix(pre.X.mean)},'Xstd',{num2cell_matrix(pre.X.std)}, ...
                  'lnKeq',{num2cell(log(r.Keq(:)))},'lnKcatf',{num2cell(log(r.Kcatf(:)))},'lnKcatr',{num2cell(log(r.Kcatr(:)))}, ...
                  'lnKM',{num2cell(full(log(r.KM(KM_i))))},'lnKA',{num2cell(full(log(r.KA(KA_i))))},'lnKI',{num2cell(full(log(r.KI(KI_i))))});
  mb_save_json_file(result,fullfile(outdir,[c.name '_matlab.json']));

  files = c.files; if ~isempty(c.constraints), files{4} = c.constraints; else, files{4} = ''; end
  manifest{end+1} = struct('name',c.name,'network',files{1},'kinetic',files{2},'state',files{3},'constraints',files{4}, ...
                           'mean',{c.mean},'std',{c.std},'use_kinetic_data',c.use_kinetic_data,'alpha',c.alpha);
end

mb_save_json_file(manifest,fullfile(outdir,'manifest.json'));


% ---------------------------------------------------


function C = num2cell_matrix(M)
% nested cell arrays (one cell per row), so that jsonencode keeps the matrix shape, also for one column
C = cell(size(M,1),1);
for i = 1:size(M,1)
  C{i} = num2cell(full(M(i,:)));
end
end
