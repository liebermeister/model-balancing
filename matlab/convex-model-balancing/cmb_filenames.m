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

filenames.model_name     = model_name;
filenames.network_sbml   = network_sbml_file;
filenames.run_dir        = [ out_DIR '/simulations/' run ];
filenames.graphics_dir   = [ out_DIR '/simulations/' run '/ps-files/' ];
filenames.data           = [ out_DIR '/simulations/' run '/data/' ];
%filenames.kinetic_true   = [ out_DIR '/simulations/' run '/data/kinetic_true.tsv'];
filenames.state_true     = [ out_DIR '/simulations/' run '/data/true'];
filenames.kinetic_data   = [ out_DIR '/simulations/' run '/data/kinetic_data.tsv'];
filenames.state_data     = [ out_DIR '/simulations/' run '/data/state_data.tsv'];
filenames.parameters_out = [ out_DIR '/simulations/' run '/data/parameters.tsv'];
filenames.states_out     = [ out_DIR '/simulations/' run '/data/'];
filenames.model_out      = [ out_DIR '/simulations/' run '/data/model.tsv'];
filenames.report         = [ out_DIR '/simulations/' run '/data/report.txt'];
filenames.options_file   = [ out_DIR '/simulations/' run '/data/options'];
filenames.options_tsv    = [ out_DIR '/simulations/' run '/data/options.tsv'];
filenames.result_file    = [ out_DIR '/simulations/' run '/data/results'];
filenames.kinetic_model    = 'kinetic_model';
filenames.metabolic_states = 'metabolic_states';

[~,~] = mkdir(filenames.graphics_dir);
[~,~] = mkdir(filenames.data);
