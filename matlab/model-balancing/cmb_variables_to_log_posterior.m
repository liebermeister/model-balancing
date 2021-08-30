function f = cmb_variables_to_log_posterior(kinetics,X,network,q_info,cmb_options,prior,data,verbose)

% f = cmb_variables_to_log_posterior(kinetics,X)
% wrapper function around  cmb_log_posterior, using the original variables (not the vector y) as inputs

clear global

pp = cmb_make_pp(network);
nn = network;
nn.kinetics = kinetics;
%modular_rate_law_haldane(nn, nn.kinetics)
q = cmb_kinetics_to_q(nn, cmb_options, q_info);
%pm(exp(q), q_info.q.names)
%nn.kinetics.Keq
%nn.kinetics.KM
%nn.kinetics.Kcatf
%nn.kinetics.Kcatr

qall = cmb_q_to_qall(q, q_info);
%pm(exp(qall), q_info.qall.names)

preposterior = cmb_prepare_posterior(prior, data, cmb_options, q_info);

[nm,ns] = size(data.X.mean);
y = cmb_qX_to_y(q,X,nm,ns);

f = cmb_log_posterior(y,pp,preposterior,data.V.mean,cmb_options,q_info,verbose);
