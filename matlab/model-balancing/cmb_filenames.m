function filenames = cmb_filenames(model_name, run, result_dir, network_sbml_file)

% filenames = cmb_filenames(model_name, result_dir, network_sbml_file)
%
% Filenames for convex estimation runs
%
% Input variables
%   model_name      (string) model name  
%   run             (string) run id
%   result_dir      Name of output directory (simple string, no directory path)    
%
% Output variables
%   filenames struct describing filenames used (constructed using [result_dir] and [run])

eval(default('result_dir','[]','network_sbml_file','[]'));

% set filenames

if ~isempty(result_dir),
  if strcmp(result_dir(end),filesep),
    result_dir = result_dir(1:end-1);
  end
  out_DIR = [result_dir filesep model_name];;
else
  out_DIR = [cmb_basedir '/../results/' model_name];
  display(sprintf('No output directory specified. I will use the directory %s',out_DIR));
end

filenames.model_name          = model_name;
filenames.network_sbml        = network_sbml_file;
filenames.run_dir             = [ out_DIR '/' run ];
filenames.graphics_dir        = [ out_DIR '/' run '/ps-files/' ];
filenames.data                = [ out_DIR '/' run '/data/' ];
filenames.model_true          = [ out_DIR '/' run '/data/model_true.tsv'];
filenames.true                = [ out_DIR '/' run '/data/true'];
filenames.kinetic_data        = [ out_DIR '/' run '/data/data_kinetic.tsv'];
filenames.state_data          = [ out_DIR '/' run '/data/data_states.tsv'];
filenames.parameters_balanced = [ out_DIR '/' run '/data/parameters_balanced.tsv'];
filenames.balanced            = [ out_DIR '/' run '/data/balanced'];
filenames.extension_model_state = 'model';
filenames.extension_state_runs  = 'states';
%filenames.balanced_model      = [ out_DIR '/' run '/data/balanced_model.tsv'];
%filenames.balanced_states     = [ out_DIR '/' run '/data/balanced_states.tsv'];
filenames.result_file         = [ out_DIR '/' run '/data/results.mat'];
filenames.report              = [ out_DIR '/' run '/data/report.txt'];
filenames.options_file        = [ out_DIR '/' run '/data/options.mat'];
filenames.options_tsv         = [ out_DIR '/' run '/data/options.tsv'];
filenames.cvxpy_json_file     = [ out_DIR '/' run '/' filesep model_name '_' run '.json'];

[~,~] = mkdir(filenames.graphics_dir);
[~,~] = mkdir(filenames.data);

