function cmb_save_results(network, data, optimal, filenames, options)

% cmb_save_results(network, data, optimal, filenames, options)
%
% save results of model balancing as SBtab files
  
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
% sbtab_document_save_to_one(sbtab_document,filenames.model_out);


% -----------------------------------------
% Save metabolic states as SBtab

clear result

for it = 1:ns,
  result.C.(data.samples{it}) = optimal.C(:,it);
  result.E.(data.samples{it}) = optimal.E(:,it);
  result.A.(data.samples{it}) = optimal.A(:,it);
  result.V.(data.samples{it}) = data.V.mean(:,it);
end

ecm_sbtab_options = struct('r', optimal.kinetics,'method','emc4cm','document_name',filenames.model_name, 'save_tolerance_ranges',0,'sbtab_attributes',struct('CalculationTime', options.calculation_time),'filename_model_state', 'kinetic_model', 'filename_state_runs','metabolic_states');

ecm_save_result_sbtab(filenames.states_out, network_optimal, result.C, result.E, result.A, ecm_sbtab_options,[],[],[],[],[],[],[],result.V);
