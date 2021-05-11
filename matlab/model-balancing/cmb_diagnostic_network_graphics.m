function cmb_diagnostic_network_graphics(network, data, optimal, cmb_options, q_info, fignum, sample, graphics_dir)

% cmb_diagnostic_network_graphics(network, data, optimal, cmb_options, q_info, fignum, sample)

eval(default('fignum','1','sample','1'));
  
kinetics_data = cmb_qall_to_kinetics(data.qall.mean,network,cmb_options,q_info);

% note that numbers in plots below refer to reaction orientations, not flux directions!

figure(fignum+1); clf
flux_signs = sign(mean(sign(data.V.mean),2)+0.5);
netgraph_concentrations(network,[], flux_signs, 1, struct('actprintnames',1));
title('Reaction orientation w.r.t flux signs')

figure(fignum); clf
subplot(2,3,1);
mm = exp(max(abs(optimal.X(:,sample)-data.X.mean(:,sample))));
netgraph_concentrations(network,optimal.X(:,sample)-data.X.mean(:,sample),[],0,struct('keep_subplot',1));
title(sprintf('met (fit/data) smpl %d f.c.<%.2f',sample,mm));

subplot(2,3,2);
mm = exp(max(abs(log(optimal.E(:,sample))-data.lnE.mean(:,sample))));
netgraph_concentrations(network,[],log(optimal.E(:,sample))-data.lnE.mean(:,sample),0,struct('keep_subplot',1));
title(sprintf('enz (fit/data) smpl %d f.c.<%.2f',sample,mm));

subplot(2,3,3);
mm = exp(max(abs(log(optimal.kinetics.Keq)-log(kinetics_data.Keq))));
netgraph_concentrations(network,[], log(optimal.kinetics.Keq)-log(kinetics_data.Keq),0,struct('keep_subplot',1));
title(sprintf('Keq (fit/data) f.c.<%.2f',mm));

subplot(2,3,4);
mm = exp(max(abs(log(optimal.kinetics.Kcatf)-log(kinetics_data.Kcatf))));
netgraph_concentrations(network,[],  log(optimal.kinetics.Kcatf)-log(kinetics_data.Kcatf),0,struct('keep_subplot',1));
title(sprintf('Kcatf (fit/data) f.c.<%.2f',mm));

subplot(2,3,5);
mm = exp(max(abs(log(optimal.kinetics.Kcatr)-log(kinetics_data.Kcatr))));
netgraph_concentrations(network,[], log(optimal.kinetics.Kcatr)-log(kinetics_data.Kcatr),0,struct('keep_subplot',1));
title(sprintf('Kcatr (fit/data) f.c.<%.2f',mm));

subplot(2,3,6);
dum = log(optimal.kinetics.KM)-log(kinetics_data.KM);
dum(isnan(dum)) = 0; 
mm = exp(max(full(abs(dum(:)))));
netgraph_concentrations(network,[], [], 0, struct('edgevalues',dum','edgestyle','fixed','keep_subplot',1));
title(sprintf('KM (fit/data) f.c.<%.2f',mm));

set(gcf,'Position', [500,100,900,800]);

if length(graphics_dir),
  print([graphics_dir 'diagnostic_errors.eps'],     sprintf('-f%d',fignum), '-depsc'); 
  print([graphics_dir 'diagnostic_flux_signs.eps'], sprintf('-f%d',fignum+1), '-depsc'); 

end