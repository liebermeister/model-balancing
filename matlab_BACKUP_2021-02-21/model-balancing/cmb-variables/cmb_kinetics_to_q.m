function q = cmb_kinetics_to_q(network, cmb_options, q_info)

[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);

r = network.kinetics;

switch cmb_options.parameterisation,
  
  case 'KV_KM_KA_KI',
    q = full([log(r.KV); log(r.KM(KM_ind)); log(r.KA(KA_ind)); log(r.KI(KI_ind))]);
  
  case 'Keq_KV_KM_KA_KI',
    ln_Keq_ind = q_info.M_qallKeq_to_qKeqind * log(r.Keq);
    q = full([ln_Keq_ind; log(r.KV); log(r.KM(KM_ind)); log(r.KA(KA_ind)); log(r.KI(KI_ind))]);
  
  otherwise error('Parameterisation not supported');

end