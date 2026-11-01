function r = mb_q_to_kinetics(q,network,mb_options,q_info)

% r = mb_q_to_kinetics(q, network, mb_options, q_info)
%
% Convert a vector q of independent log kinetic constants into a kinetics struct
% (format as in network.kinetics, Metabolic Network Toolbox); kcat values are computed from KV and Keq

[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);

r = network.kinetics;

r.Keq = exp(q_info.M_qKeqind_to_qallKeq * q(q_info.q.index.Keq));

r.KV         = exp(q(q_info.q.index.KV));
r.KM(KM_ind) = exp(q(q_info.q.index.KM));
r.KA(KA_ind) = exp(q(q_info.q.index.KA));
r.KI(KI_ind) = exp(q(q_info.q.index.KI));

[r.Kcatf, r.Kcatr] = modular_KV_Keq_to_kcat(network.N,r);
