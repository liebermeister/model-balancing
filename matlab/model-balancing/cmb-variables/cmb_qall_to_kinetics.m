function r = cmb_qall_to_kinetics(qall,network,cmb_options,q_info)

[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);

r = network.kinetics;

switch cmb_options.parameterisation,
  case 'Keq_KV_KM_KA_KI',
    r.Keq = exp(q_info.M_qKeqind_to_qallKeq * qall(q_info.qall.index.Keq));
end

r.KV         = exp(qall(q_info.qall.index.KV));
r.KM(KM_ind) = exp(qall(q_info.qall.index.KM));
r.KA(KA_ind) = exp(qall(q_info.qall.index.KA));
r.KI(KI_ind) = exp(qall(q_info.qall.index.KI));
r.Kcatf      = exp(qall(q_info.qall.index.Kcatf));
r.Kcatr      = exp(qall(q_info.qall.index.Kcatr));
