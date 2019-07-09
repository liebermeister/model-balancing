function [network, q_info, data, c_init, kinetic_data, state_data] = e_coli_davidi_2016_model_and_data(cmb_options)

% --------------------------------------------------------------
% Set model name and input filenames

[model_file, kinetic_data_file] = e_coli_glucose_kinetic_data_filenames(cmb_options);
state_data_files                = e_coli_davidi_2016_data_filenames(cmb_options);

replace_ids_in_network = {{'R04779','R00756'}}; % compatibility with davidi annotations

[network, kinetic_data, state_data] = cmb_load_model_and_data(model_file, kinetic_data_file, state_data_files,replace_ids_in_network);

% Flux unit conversion from mmol/gdw/s to mM/s (use 1gdw = 4 g wet weight = 4 ml)
% 1 mmol/[1 gdw * 1 s] = mmol / [0.004 l * 1 s] = 252 mM/s

state_data.flux_data.Mean = 252 * state_data.flux_data.Mean;
state_data.flux_data.Std  = 252 * state_data.flux_data.Std;
state_data.flux_data.Unit = 'mM/s';

% since no standard deviations were given ..
state_data.flux_data.Std = 0.1 * abs(state_data.flux_data.Mean);
state_data.flux_data.Std(state_data.flux_data.Std==0) = 0.01 * max(abs(state_data.flux_data.Mean(:)));

% Metabolite unit conversion from uM to mM

state_data.metabolite_data.Mean = 0.001 * state_data.metabolite_data.Mean;
state_data.metabolite_data.Std  = 0.001 * state_data.metabolite_data.Std;
state_data.metabolite_data.Unit = 'mM';

% Enzyme unit conversion from mmol/gCDW to mM
% 1 mmol/[1 gdw] = mmol / [0.004 l] = 252 mM

state_data.enzyme_data.Mean = 252 * state_data.enzyme_data.Mean;
state_data.enzyme_data.Std  = 252 * state_data.enzyme_data.Std;
state_data.enzyme_data.Unit = 'mM';

data = cmb_state_data_to_data(state_data, cmb_options);

% Complete flux data by using projection (MOVE THIS ELSEWHERE???)

[nr,ns] = size(state_data.flux_data.Mean);
v_sign = ones(nr,1);
V     = state_data.flux_data.Mean;
V_std = state_data.flux_data.Std;
for it = 1:ns,
  [V_proj(:,it), ~, V_proj_std(:,it)] = project_fluxes(network.N, find(network.external), V(:,it), V_std(:,it),v_sign);
end

data.V.mean = V_proj;
data.V.std  = V_proj_std;

% Initial parameters are a bit stupid... 

[nr,nm,nx,KM_indices,KA_indices,KI_indices,nKM,nKA,nKI] = network_numbers(network);

q_info = cmb_define_parameterisation(network, cmb_options); 

nn             = network;
nn.kinetics    = data_integration_median_values_to_kinetics(network.kinetics, kinetic_data);
data.qall.mean = nan * ones(q_info.qall.number,1);
data.qall.mean(q_info.qall.index.Keq)   = log(nn.kinetics.Keq);
data.qall.mean(q_info.qall.index.KM)    = log(nn.kinetics.KM(KM_indices));
data.qall.mean(q_info.qall.index.KA)    = log(nn.kinetics.KA(KA_indices));
data.qall.mean(q_info.qall.index.KI)    = log(nn.kinetics.KI(KI_indices));
data.qall.mean(q_info.qall.index.Kcatf) = log(nn.kinetics.Kcatf);
data.qall.mean(q_info.qall.index.Kcatr) = log(nn.kinetics.Kcatr);
data.qall.std  = log(cmb_options.data_kin_geom_std) * ones(q_info.qall.number,1);
data.qall.std(isnan(data.qall.mean)) = nan;

% Initial state for optimisation (using parameters from network.kinetics and concentrations = 1)

c_init = ones(nm,1);
