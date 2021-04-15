function filenames = cmb_filenames(model_name, prb, run, result_dir, network_sbml_file, flag_artificial)

% filenames = cmb_filenames(model_name, prb, run, result_dir, network_sbml_file, flag_artificial)
%
% Filenames for convex estimation runs
%
% Input variables
%   model_name      (string) model name  
%   prb             (string) problem id (should not refer to alpha value or other optimisation parameters)
%   run             (string) run id (like sim, but may also refer to alpha value or other optimisation parameters)
%   result_dir      Name of output directory (simple string, no directory path)    
%
% Output variables
%   filenames struct describing filenames used (constructed using [result_dir] and [run])

eval(default('result_dir','[]','network_sbml_file','[]','flag_artificial','0'));

% set filenames

if ~isempty(result_dir),
  if strcmp(result_dir(end),filesep),
    result_dir = result_dir(1:end-1);
  end
  out_DIR = [result_dir filesep model_name];;
else
  out_DIR = [cmb_basedir '/results/' model_name];
  display(sprintf('No output directory specified. I will use the directory %s',out_DIR));
end

filenames.model_name            = model_name;
filenames.network_sbml          = network_sbml_file;
run_dir                         = [ out_DIR '/' run ];
data_dir                        = [ out_DIR '/' run '/data/' ];
filenames.graphics_dir          = [ out_DIR '/' run '/ps-files/' ];
if flag_artificial,
  %filenames.model_true           = [ out_DIR '/' run '/data/model_true.tsv'];
  filenames.true                  = [ out_DIR '/' run '/data/true'];
end
%filenames.parameters_balanced  = [ out_DIR '/' run '/data/balanced_parameters.tsv'];
filenames.balanced              = [ out_DIR '/' run '/data/balanced'];
filenames.extension_model_state = 'model';
filenames.extension_state_runs  = 'states';
filenames.kinetic_data          = [ out_DIR '/' run '/data/data_kinetic.tsv'];
filenames.state_data            = [ out_DIR '/' run '/data/data_states.tsv'];
%filenames.balanced_model       = [ out_DIR '/' run '/data/balanced_model.tsv'];
%filenames.balanced_states      = [ out_DIR '/' run '/data/balanced_states.tsv'];
filenames.results_mat           = [ out_DIR '/' run '/data/results.mat'];
filenames.report_txt            = [ out_DIR '/' run '/data/report.txt'];
%filenames.options_file         = [ out_DIR '/' run '/data/options.mat'];
filenames.options_sbtab         = [ out_DIR '/' run '/data/options.tsv'];
filenames.model_and_data_json   = [ out_DIR '/' run '/' filesep model_name '_' prb '.json'];

[~,~] = mkdir(run_dir);
[~,~] = mkdir(data_dir);
[~,~] = mkdir(filenames.graphics_dir);
