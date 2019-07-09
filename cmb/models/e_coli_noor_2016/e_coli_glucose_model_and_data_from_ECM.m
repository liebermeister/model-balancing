function [network, q_info, data, c_init] = e_coli_model_and_data_from_ECM_result(cmb_options)

% --------------------------------------------------------------
% Load E. coli model and input data from ECM paper

DIR = [cmb_basedir '/resources/escherichia_coli/noor_2016'];
load([DIR '/ecoli_ccm_ProteinComposition_Haverkorn_model_data']);
load([DIR '/ecoli_ccm_ProteinComposition_Haverkorn_options']);
load([DIR '/ecoli_ccm_ProteinComposition_Haverkorn']);


% --------------------------------------------------------------
% Data

[nr,nm] = network_numbers(network);

q_info  = cmb_define_parameterisation(network, cmb_options); 

cmb_options.ns = 1;

for it = 1:cmb_options.ns,
  data.samples{it,1} = ['S' num2str(it)];
end

data.X.mean = log(c_data);
data.X.std  = log(cmb_options.data_C_geom_std) * ones(nm,1);
data.E.mean = u_data;
data.E.std  = [cmb_options.data_E_geom_std-1] * data.E.mean;
data.V.mean = v;
data.V.std  = v_std;

% use kinetic constants in structure "network" (loaded from ECM result)
q_all_from_ECM = cmb_q_to_qall(cmb_kinetics_to_q(network, cmb_options, q_info),q_info);

display('Using kinetic constants from E. coli ECM model as surrogate kinetic data');
data.qall.mean = q_all_from_ECM;
data.qall.std  = log(cmb_options.data_kin_geom_std) * ones(q_info.qall.number,1);

% Initial state for optimisation (from ECM)

c_init = c.emc4cm;

