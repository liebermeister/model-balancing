function r = cmb_q_to_kinetics(q,network,cmb_options,q_info)

[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);

r = network.kinetics;

switch cmb_options.parameterisation,
  case 'Keq_KV_KM_KA_KI',
    r.Keq = exp(q_info.M_qKeqind_to_qallKeq * q(q_info.q.index.Keq) );
end

r.KV         = exp(q(q_info.q.index.KV));
r.KM(KM_ind) = exp(q(q_info.q.index.KM));
r.KA(KA_ind) = exp(q(q_info.q.index.KA));
r.KI(KI_ind) = exp(q(q_info.q.index.KI));

[r.Kcatf, r.Kcatr] = modular_KV_Keq_to_kcat(network.N,r);
