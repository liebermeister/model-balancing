function [network, q_info, data, c_init, kinetic_data, state_data] = e_coli_glucose_model_and_data(cmb_options)

%% Generic function for reading models and data - does not yield a good initial state! 

% --------------------------------------------------------------
% Set model name and input filenames

[model_file, kinetic_data_file] = e_coli_glucose_kinetic_data_filenames(cmb_options);
state_data_files                = e_coli_glucose_data_filenames(cmb_options);

[network, kinetic_data, state_data] = cmb_load_model_and_data(model_file, kinetic_data_file, state_data_files);

% Fluces from Chen et al (2010), Fig. 1 (some rows were duplicated) 
% Flux unit:   (mmol/g dry weight/h)
% Growth rate: Aerobic:   0.58 7 0.01 (1/h) 
%
% Flux unit conversion from mmol/gdw/h to mM/s (use 1gdw = 4 g wet weight = 4 ml)
% 1 mmol/[1 gdw * 1 h] = mmol / [0.004 l * 3600 s] = 0.07 mM/s

state_data.flux_data.Mean = 0.07 * state_data.flux_data.Mean;
state_data.flux_data.Std  = 0.07 * state_data.flux_data.Std;
state_data.flux_data.Unit = 'mM/s';

data = cmb_state_data_to_data(state_data, cmb_options);

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

% use data value for concentrations and complete missing values
c_init = state_data.metabolite_data.Mean;
c_init(isnan(c_init)) = 0.01;
