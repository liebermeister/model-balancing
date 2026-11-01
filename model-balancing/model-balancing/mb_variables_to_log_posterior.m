function f = mb_variables_to_log_posterior(kinetics,X,network,q_info,mb_options,prior,data,verbose)

% f = mb_variables_to_log_posterior(kinetics,X,network,q_info,mb_options,prior,data,verbose)
%
% Wrapper function around mb_log_posterior, using the original variables (not the vector y) as inputs

clear global

pp = mb_make_pp(network);
nn = network;
nn.kinetics = kinetics;
%modular_rate_law_haldane(nn, nn.kinetics)
q = mb_kinetics_to_q(nn, mb_options, q_info);
%pm(exp(q), q_info.q.names)
%nn.kinetics.Keq
%nn.kinetics.KM
%nn.kinetics.Kcatf
%nn.kinetics.Kcatr

qall = mb_q_to_qall(q, q_info);
%pm(exp(qall), q_info.qall.names)

preposterior = mb_prepare_posterior(prior, data, mb_options, q_info);

[nm,ns] = size(data.X.mean);
y = mb_qX_to_y(q,X,nm,ns);

f = mb_log_posterior(y,pp,preposterior,data.V.mean,mb_options,q_info,verbose);
