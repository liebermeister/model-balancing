function res = cmb_qX_to_variables(q,X,V,network,cmb_options,q_info,ecm_score,pp)

[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);

ns = size(X,2);

kinetics = cmb_q_to_kinetics(q,network,cmb_options,q_info);

pp.network.kinetics = kinetics;

for it = 1:ns,
  pp.v = V(:,it);
  [~, E(:,it)] = ecm_get_score(ecm_score,X(:,it),pp);
end

res.kinetics = kinetics;
res.q        = q;
res.X        = X;
res.C        = exp(X);
res.E        = E;
res.E(find(V==0)) = 0;

% reaction affinity in kJ/mol, with signs
for it = 1:ns,
  res.A(:,it) = RT * [log(kinetics.Keq) - network.N' * X(:,it)];
end

% check whether solutions yield the right fluxes
%  network.kinetics = kinetics;
%  for it = 1:ns,
%    network.kinetics.u = res.E(:,it);
%    [V(:,it), network_velocities(res.C(:,it),network)]
%  end
