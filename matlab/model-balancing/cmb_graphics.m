function cmb_graphics(prior, data, optimal, true, cmb_options, q_info, graphics_dir, kapp_max)

% ---------------------
% set some variables
  
eval(default('true','[]'));

qall_names   = q_info.qall.names;
qall_data    = data.qall.mean;
qall_optimal = cmb_q_to_qall(optimal.q,q_info);
index        = q_info.qall.index;
nr           = length(index.KV);

% ---------------------
% Graphics: plot data vs fitted

display(sprintf('Writing graphics files with basename %s',graphics_dir));

filenames.x        = [graphics_dir 'data_x.eps'];
filenames.e        = [graphics_dir 'data_e.eps'];
filenames.KV       = [graphics_dir 'data_KV.eps'];
filenames.KM       = [graphics_dir 'data_KM.eps'];
filenames.Keq      = [graphics_dir 'data_Keq.eps'];
filenames.Kcatf    = [graphics_dir 'data_Kcatf.eps'];
filenames.Kcatr    = [graphics_dir 'data_Kcatr.eps'];
filenames.Kappmaxf = [graphics_dir 'data_Kappmaxf.eps'];
filenames.Kappmaxr = [graphics_dir 'data_Kappmaxr.eps'];

fit_color  = [1 0 0];
pred_color = [1 0 1];
met_color  = [0 0 0.8];
enz_color  = [.6 0.2 0.1];
kapp_color = [.7 0 .7];

switch cmb_options.use_kinetic_data,
  case 'all', 
    keq_color = fit_color; kcat_color = fit_color; km_color = fit_color;
  case 'only_Keq_data',
    keq_color = fit_color; kcat_color = pred_color; km_color = pred_color;
  case 'none',
    keq_color = pred_color; kcat_color = pred_color; km_color = pred_color;
end

figure(1); 
scatter_plot(data.X.mean,optimal.X,[],met_color,'Metabolite levels [mM] (data)', 'Metabolite levels [mM] (fit)');

figure(2); 
%scatter_plot(log(data.E.mean),log(optimal.E),[],enz_color,'Enzyme levels [mM] (data)', 'Enzyme levels [mM] (fit)');
scatter_plot(data.lnE.mean,log(optimal.E),[],enz_color,'Enzyme levels [mM] (data)', 'Enzyme levels [mM] (fit)');

print([filenames.x], '-f1', '-depsc'); 
print([filenames.e], '-f2', '-depsc'); 

figure(3); 
scatter_plot(qall_data(index.KV),qall_optimal(index.KV),qall_names(index.KV),fit_color,'K_{V} values [1/s] (data)', 'K_{V} values [1/s] (fit)');
    
figure(4);
% % mark KM value if the minimum concentration, over all samples, is at least twice as large as KM value
% % (consider optimal values for this criterion)
% % log metabolite concentrations, minimum across all samples
% c_min_over_samples = min(exp(optimal.X),[],2);
% C_min_mat = repmat(c_min_over_samples',nr,1);
% mark_indices = find(C_min_mat(q_info.KM_matrix_indices) > 2 * exp(qall_optimal(index.KM)));

scatter_plot(qall_data(index.KM),qall_optimal(index.KM),qall_names(index.KM),km_color,'K_{M} values [mM] (data)', 'K_{M} values [mM] (fit)');

print([filenames.KV], '-f3', '-depsc'); 
print([filenames.KM], '-f4', '-depsc'); 

switch cmb_options.parameterisation,
  case 'Keq_KV_KM_KA_KI',
    
    figure(5); 
    scatter_plot(qall_data(index.Keq),qall_optimal(index.Keq),qall_names(index.Keq),keq_color,'K_{eq} values [unitless] (data)', 'K_{eq} values [unitless] (fit)');
      
    figure(6); 
    scatter_plot(qall_data(index.Kcatf),qall_optimal(index.Kcatf),qall_names(index.Kcatf),kcat_color,'k_{cat}^{+} values [1/s] (data)', 'k_{cat}^{+} values [1/s] (fit)');
    
    figure(7); 
    scatter_plot(qall_data(index.Kcatr),qall_optimal(index.Kcatr),qall_names(index.Kcatr),kcat_color,'k_{cat}^{-} values [1/s] (data)', 'k_{cat}^{-} values [1/s] (fit)');
      
    figure(8); 
    scatter_plot(qall_data(index.Kcatf),log(kapp_max.forward),qall_names(index.Kcatf),kapp_color,'k_{cat}^{+} values [1/s] (data)', 'k_{app,max}^{+} values [1/s] (fit)');

    
    print([filenames.Keq],   '-f5', '-depsc'); 
    print([filenames.Kcatf], '-f6', '-depsc'); 
    print([filenames.Kcatr], '-f7', '-depsc'); 
    print([filenames.Kappmaxf], '-f8', '-depsc'); 

    %% only if some reverse kapp values have been determined:
    if sum(isfinite(kapp_max.reverse)),
      figure(9);
      scatter_plot(qall_data(index.Kcatr),log(kapp_max.reverse),qall_names(index.Kcatr),kapp_color,'k_{cat}^{-} values [1/s] (data)', 'k_{app,max}^{-} values [1/s] (fit)');
      print([filenames.Kappmaxr], '-f9', '-depsc'); 
    end

end
  
% -----------------------------------------------------------
% Graphics: plot true vs fitted (only in the case of artificial data)

if length(true),

  filenames.x        = [graphics_dir 'true_x.eps'];
  filenames.e        = [graphics_dir 'true_e.eps'];
  filenames.KV       = [graphics_dir 'true_KV.eps'];
  filenames.KM       = [graphics_dir 'true_KM.eps'];
  filenames.Keq      = [graphics_dir 'true_Keq.eps'];
  filenames.Kcatf    = [graphics_dir 'true_Kcatf.eps'];
  filenames.Kcatr    = [graphics_dir 'true_Kcatr.eps'];
  filenames.Kappmaxf = [graphics_dir 'true_Kappmaxf.eps'];
  filenames.Kappmaxr = [graphics_dir 'true_Kappmaxr.eps'];

  qall_true = cmb_q_to_qall(true.q,q_info);

  figure(11); 
  scatter_plot(true.X,optimal.X,[],met_color,'Metabolite levels [mM] (true)','Metabolite levels [mM] (fit)');
  
  figure(12); 
  scatter_plot(log(true.E),log(optimal.E),[],enz_color,'Enzyme levels [mM] (true)','Enzyme levels [mM] (fit)');

  figure(13); 
  scatter_plot(qall_true(index.KV),qall_optimal(index.KV),qall_names(index.KV),fit_color,'K_{V} values [1/s] (true)', 'K_{V} values [1/s] (fit)');
  
  figure(14); 
  % % mark KM value if the minimum concentration, over all samples, is at least twice as large as KM value
  % % (consider optimal values for this criterion)
  % % log metabolite concentrations, minimum across all samples
  % c_min_over_samples = min(exp(optimal.X),[],2);
  % C_min_mat = repmat(c_min_over_samples',nr,1);
  % mark_indices = find(C_min_mat(q_info.KM_matrix_indices) > 2 * exp(qall_optimal(index.KM)));

  scatter_plot(qall_true(index.KM),qall_optimal(index.KM),qall_names(index.KM),km_color,'K_{M} values [mM] (true)', 'K_{M} values [mM] (fit)');
  
  print([filenames.x], '-f11', '-depsc'); 
  print([filenames.e], '-f12', '-depsc'); 
  print([filenames.KV], '-f13', '-depsc'); 
  print([filenames.KM], '-f14', '-depsc'); 
    

  switch cmb_options.parameterisation,
    case 'Keq_KV_KM_KA_KI',
      
      figure(15); 
      scatter_plot(qall_true(index.Keq),qall_optimal(index.Keq),qall_names(index.Keq),keq_color,'K_{eq} values [unitless] (true)', 'K_{eq} values [unitless] (fit)');

      figure(16); 
      scatter_plot(qall_true(index.Kcatf),qall_optimal(index.Kcatf),qall_names(index.Kcatf),kcat_color,'k_{cat}^{+} values [1/s] (true)', 'k_{cat}^{+} values [1/s] (fit)');
      
      figure(17); 
      scatter_plot(qall_true(index.Kcatr),qall_optimal(index.Kcatr),qall_names(index.Kcatr),kcat_color,'k_{cat}^{-} values [1/s] (true)', 'k_{cat}^{-} values [1/s] (fit)');

      figure(18); 
      scatter_plot(qall_true(index.Kcatf),log(kapp_max.forward),qall_names(index.Kcatf),kapp_color,'k_{cat}^{+} values [1/s] (true)', 'k_{app,max}^{+} values [1/s] (fit)');

      print([filenames.Keq],   '-f15', '-depsc'); 
      print([filenames.Kcatf], '-f16', '-depsc'); 
      print([filenames.Kcatr], '-f17', '-depsc'); 
      print([filenames.Kappmaxf], '-f18', '-depsc'); 

      %% only if some reverse kapp values have been determined
      if sum(isfinite(kapp_max.reverse)),
        figure(19); 
        scatter_plot(qall_true(index.Kcatr),log(kapp_max.reverse),qall_names(index.Kcatr),kapp_color,'k_{cat}^{-} values [1/s] (true)', 'k_{app,max}^{-} values [1/s] (fit)');
      print([filenames.Kappmaxr], '-f19', '-depsc'); 
      end
      
  end

end

% -----------------------------------------------------------
% Graphics: plot true vs data (only in the case of artificial data)

if cmb_options.plot_true_vs_data,
if length(true),

  filenames.x        = [graphics_dir 'trueVSdata_x.eps'];
  filenames.e        = [graphics_dir 'trueVSdata_e.eps'];
  filenames.KV       = [graphics_dir 'trueVSdata_KV.eps'];
  filenames.KM       = [graphics_dir 'trueVSdata_KM.eps'];
  filenames.Keq      = [graphics_dir 'trueVSdata_Keq.eps'];
  filenames.Kcatf    = [graphics_dir 'trueVSdata_Kcatf.eps'];
  filenames.Kcatr    = [graphics_dir 'trueVSdata_Kcatr.eps'];

  switch cmb_options.use_kinetic_data
    case 'all',
      kin_color = fit_color;
      keq_color = fit_color;
    case 'only_Keq_data'
      kin_color = pred_color;
      keq_color = fit_color;
    case 'none'
      kin_color = pred_color;
      keq_color = pred_color;
  end
  
  qall_true = cmb_q_to_qall(true.q,q_info);

  figure(21); 
  scatter_plot(true.X,data.X.mean,[],met_color,'Metabolite levels [mM] (true)','Metabolite levels [mM] (data)');
  
  figure(22); 
  scatter_plot(log(true.E),log(data.E.mean),[],enz_color,'Enzyme levels [mM] (true)','Enzyme levels [mM] (data)');

  figure(23); 
  scatter_plot(qall_true(index.KV),qall_data(index.KV),qall_names(index.KV),kin_color,'K_{V} values [1/s] (true)', 'K_{V} values [1/s] (data)');
  
  figure(24); 
  % % mark KM value if the minimum concentration, over all samples, is at least twice as large as KM value
  % % (consider true values for this criterion)
  % % log metabolite concentrations, minimum across all samples
  % c_min_over_samples = min(exp(true.X),[],2);
  % C_min_mat = repmat(c_min_over_samples',nr,1);
  % mark_indices = find(C_min_mat(q_info.KM_matrix_indices) > 2 * exp(qall_true(index.KM)));

  scatter_plot(qall_true(index.KM),qall_data(index.KM),qall_names(index.KM),kin_color,'K_{M} values [mM] (true)', 'K_{M} values [mM] (data)');
  
  print([filenames.x], '-f21', '-depsc'); 
  print([filenames.e], '-f22', '-depsc'); 
  print([filenames.KV], '-f23', '-depsc'); 
  print([filenames.KM], '-f24', '-depsc'); 
    

  switch cmb_options.parameterisation,
    case 'Keq_KV_KM_KA_KI',
      
      figure(25); 
      scatter_plot(qall_true(index.Keq),qall_data(index.Keq),qall_names(index.Keq),keq_color,'K_{eq} values [unitless] (true)', 'K_{eq} values [unitless] (data)');

      figure(26); 
      scatter_plot(qall_true(index.Kcatf),qall_data(index.Kcatf),qall_names(index.Kcatf),kin_color,'k_{cat}^{+} values [1/s] (true)', 'k_{cat}^{+} values [1/s] (data)');
      
      figure(27); 
      scatter_plot(qall_true(index.Kcatr),qall_data(index.Kcatr),qall_names(index.Kcatr),kin_color,'k_{cat}^{-} values [1/s] (true)', 'k_{cat}^{-} values [1/s] (data)');

      print([filenames.Keq],   '-f25', '-depsc'); 
      print([filenames.Kcatf], '-f26', '-depsc'); 
      print([filenames.Kcatr], '-f27', '-depsc'); 

  end

end
end

% ================================================

function scatter_plot(xvalues,yvalues,names,color,xtitle,ytitle,mark_indices)

eval(default('mark_indices','[]'));

ind_ok = find(isfinite(xvalues).*isfinite(yvalues));
xvalues = xvalues(ind_ok); 
yvalues = yvalues(ind_ok); 

% all data are assumed to be on natural log scale  

geom_deviation = exp(sqrt(nanmean([yvalues(:) - xvalues(:)].^2)));

ind = find(isfinite(yvalues) .* isfinite(xvalues));
if length(ind)>1,
  linear_correlation = corrcoef(yvalues(ind),xvalues(ind)); 
  linear_correlation = linear_correlation(1,2);
else 
  linear_correlation = nan;
end

mmin = min([exp(xvalues);exp(yvalues)]);
mmax = max([exp(xvalues);exp(yvalues)]);

mmin = 10.^floor(log10(mmin));
mmax = 10.^ceil(log10(mmax));

clf;
plot([mmin,mmax],[mmin,mmax],'-','Color',[.5 .5 .5]); hold on
plot(exp(xvalues),exp(yvalues),'.','Markersize',20,'Color',color);
plot(exp(xvalues(mark_indices)),exp(yvalues(mark_indices)),'o','Markersize',10,'Color',color);
%if length(names), text(exp(xvalues),exp(yvalues),names); end

axis equal; axis square; axis tight
set(gca,'XScale','log','YScale','log'); set(gca,'fontsize',16);
title(sprintf('GeomDev: %2.2f  CorrCoeff: %2.2f',geom_deviation,linear_correlation),'FontSize',12);
xlabel(xtitle); ylabel(ytitle);

a = axis; 
XTick = 10.^[floor(log10(min(axis))):ceil(log10(max(axis)))];
XTickLabels = cellstr(num2str(round(log10(XTick(:))), '10^{%d}'));
set(gca,'xtick',XTick,'ytick',XTick,'xticklabels',XTickLabels,'yticklabels',XTickLabels);