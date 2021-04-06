function cmb_save_results(network, data, bounds, optimal, filenames, cmb_options, options, true)

% cmb_save_results(network, data, optimal, filenames, cmb_options, options)
%
% save results of model balancing as SBtab files
% 
% filenames: struct with fields:
%   .model_name            string (used in SBtab as "document_name")
%   .parameters_out        SBtab balanced kinetic parameters
%   .true                  SBtab model and data (true)     - basename
%   .balanced              SBtab model and data (balanced) - basename
%   .extension_model_state SBtab model - extension for model file
%   .extension_state_runs  SBtab model - extension for state data file
%   .options_sbtab         options output file (SBtab) 

eval(default('true','[]','options','struct'));

ns = size(optimal.C,2);
  
network_optimal = network;
network_optimal.kinetics = optimal.kinetics;

sbtab_options = struct('use_sbml_ids',0,'verbose',0,'modular_rate_law_kinetics', 0, 'write_concentrations',0,'write_enzyme_concentrations',0,'document_name',filenames.model_name,'flag_minimal_output',0,'value_column_name', 'Value');


% -----------------------------------------
% Save kinetic parameters as SBtab

% Omitted, because redundant with files written below
% quantity_table = modular_rate_law_to_sbtab(network_optimal,[],sbtab_options);
% sbtab_table_save(quantity_table,struct('filename',filenames.parameters_out,'verbose',1));


% -----------------------------------------
% Save parameterised model as SBtab

% Omitted, because redundant with files written below
% sbtab_document = network_to_sbtab(network_optimal, sbtab_options);
% sbtab_document_save_to_one(sbtab_document,[filenames.balanced '_' filenames.extension_model_state]);


% -----------------------------------------
% Save metabolic states as SBtab

clear result

for it = 1:ns,
  result.C.(data.samples{it}) = optimal.C(:,it);
  result.E.(data.samples{it}) = optimal.E(:,it);
  result.DeltaG.(data.samples{it}) = -optimal.A(:,it);
  result.V.(data.samples{it}) = data.V.mean(:,it);
end

ecm_sbtab_options = struct('r', optimal.kinetics,'method','emc4cm','document_name',filenames.model_name, 'save_tolerance_ranges',0,'sbtab_attributes',struct('DocumentName', 'CMB result', 'Document', 'CMBresult', 'RelaxationAlpha', sprintf('%f',cmb_options.enzyme_score_alpha), 'CalculationTime', sprintf('%f s',options.calculation_time)),'filename_model_state', filenames.extension_model_state, 'filename_state_runs',filenames.extension_state_runs);

ecm_save_result_sbtab(filenames.balanced, network_optimal, result.C, result.E, result.DeltaG, ecm_sbtab_options,bounds.conc_min, bounds.conc_max,[],[],[],[],[],result.V);


% -----------------------------------------
% For artificial data: Save "true" metabolic states as SBtab

if length(true),

  network_true = network;
  network_true.kinetics = true.kinetics;

  clear result
  
  for it = 1:ns,
    result.C.(data.samples{it}) = exp(true.X(:,it));
    result.E.(data.samples{it}) = true.E(:,it);
    result.A.(data.samples{it}) = true.A_forward(:,it) .* sign(true.V(:,it));
    result.V.(data.samples{it}) = true.V(:,it);
  end
  
  ecm_sbtab_options = struct('r', optimal.kinetics,'method','emc4cm','document_name',filenames.model_name, 'save_tolerance_ranges',0,'sbtab_attributes',struct('DocumentName', 'CMB result', 'Document', 'CMBresult', 'CalculationTime', sprintf('%f s',options.calculation_time)),'filename_model_state', filenames.extension_model_state, 'filename_state_runs', filenames.extension_state_runs);

  ecm_save_result_sbtab(filenames.true, network_true, result.C, result.E, result.A, ecm_sbtab_options, bounds.conc_min, bounds.conc_max,[],[],[],[],[],result.V);

end

% ----------------------------------------
% save cmb_options as sbtab file

cmb_options_sbtab = options_to_sbtab(cmb_options,struct('TableName','Options for convex model balancing','TableID','ConfigureModelBalancing','Method','convex-model-balancing'));

if length(filenames.options_sbtab)
  sbtab_table_save(cmb_options_sbtab,struct('filename',filenames.options_sbtab));
end