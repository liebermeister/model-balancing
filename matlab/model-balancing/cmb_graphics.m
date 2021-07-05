function [res, res_true, res_true_data] = cmb_graphics(prior, data, optimal, true, cmb_options, q_info, graphics_dir, kapp_max, show_graphics)

eval(default('graphics_dir','[]','show_graphics','1'));

res = []; res_true = []; res_true_data = [];

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

if length(graphics_dir),
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
end

fit_color  = [1 0 0];
pred_color = [1 0 1];
met_color  = [0 0 0.8];
enz_color  = [.6 0.2 0.1];
kapp_color = [.7 0 .7];

switch cmb_options.use_kinetic_data,
  case 'all', 
    keq_color = fit_color; kcat_color = fit_color; km_color = fit_color;
  case 'only_Keq',
    keq_color = fit_color; kcat_color = pred_color; km_color = pred_color;
  case 'none',
    keq_color = pred_color; kcat_color = pred_color; km_color = pred_color;
end

% figure(1); 
res.c = scatter_plot(data.X.mean,optimal.X,[],met_color,'Metabolite levels [mM] (data)', 'Metabolite levels [mM] (fit)', [], show_graphics,1,filenames.x);

% figure(2); 
%res = scatter_plot(log(data.E.mean),log(optimal.E),[],enz_color,'Enzyme levels [mM] (data)', 'Enzyme levels [mM] (fit)', [], show_graphics,2);
res.e = scatter_plot(data.lnE.mean,log(optimal.E),[],enz_color,'Enzyme levels [mM] (data)', 'Enzyme levels [mM] (fit)', [], show_graphics,2,filenames.e);

% figure(3); 
res.KV = scatter_plot(qall_data(index.KV),qall_optimal(index.KV),qall_names(index.KV),fit_color,'K_{V} values [1/s] (data)', 'K_{V} values [1/s] (fit)', [], show_graphics,3,filenames.KV);
    
% figure(4);
% % mark KM value if the minimum concentration, over all samples, is at least twice as large as KM value
% % (consider optimal values for this criterion)
% % log metabolite concentrations, minimum across all samples
% c_min_over_samples = min(exp(optimal.X),[],2);
% C_min_mat = repmat(c_min_over_samples',nr,1);
% mark_indices = find(C_min_mat(q_info.KM_matrix_indices) > 2 * exp(qall_optimal(index.KM)));

res.KM = scatter_plot(qall_data(index.KM),qall_optimal(index.KM),qall_names(index.KM),km_color,'K_{M} values [mM] (data)', 'K_{M} values [mM] (fit)', [], show_graphics,4,filenames.KM);

switch cmb_options.parameterisation,
  case 'Keq_KV_KM_KA_KI',
    
    % figure(5); 
    res.Keq = scatter_plot(qall_data(index.Keq),qall_optimal(index.Keq),qall_names(index.Keq),keq_color,'K_{eq} values [unitless] (data)', 'K_{eq} values [unitless] (fit)', [], show_graphics,5,filenames.Keq);
      
    % figure(6); 
    res.Kcatf = scatter_plot(qall_data(index.Kcatf),qall_optimal(index.Kcatf),qall_names(index.Kcatf),kcat_color,'k_{cat}^{+} values [1/s] (data)', 'k_{cat}^{+} values [1/s] (fit)', [], show_graphics,6,filenames.Kcatf);
    
    % figure(7); 
    res.Kcatr = scatter_plot(qall_data(index.Kcatr),qall_optimal(index.Kcatr),qall_names(index.Kcatr),kcat_color,'k_{cat}^{-} values [1/s] (data)', 'k_{cat}^{-} values [1/s] (fit)', [], show_graphics,7,filenames.Kcatr);
      
    % figure(8); 
    res.Kappf = scatter_plot(qall_data(index.Kcatf),log(kapp_max.forward),qall_names(index.Kcatf),kapp_color,'k_{cat}^{+} values [1/s] (data)', 'k_{app,max}^{+} values [1/s] (fit)', [], show_graphics,8,filenames.Kappmaxf);

    % figure(9);
    %% plot only if some reverse kapp values have been determined
    if sum(isfinite(kapp_max.reverse)),
      res.Kappr = scatter_plot(qall_data(index.Kcatr),log(kapp_max.reverse),qall_names(index.Kcatr),kapp_color,'k_{cat}^{-} values [1/s] (data)', 'k_{app,max}^{-} values [1/s] (fit)', [], show_graphics,9,filenames.Kappmaxr);
    elseif show_graphics,
      figure(9);
      title('No reverse Kapp values determined');
      print(filenames.Kappmaxr, '-f9', '-depsc'); 
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

  % figure(11); 
  res_true.c = scatter_plot(true.X,optimal.X,[],met_color,'Metabolite levels [mM] (true)','Metabolite levels [mM] (fit)', [], show_graphics,11,filenames.x);
  
  % figure(12); 
  res_true.e = scatter_plot(log(true.E),log(optimal.E),[],enz_color,'Enzyme levels [mM] (true)','Enzyme levels [mM] (fit)', [], show_graphics,12,filenames.e);

  % figure(13); 
  res_true.KV = scatter_plot(qall_true(index.KV),qall_optimal(index.KV),qall_names(index.KV),fit_color,'K_{V} values [1/s] (true)', 'K_{V} values [1/s] (fit)', [], show_graphics,13,filenames.KV);
  
  % figure(14); 
  % % mark KM value if the minimum concentration, over all samples, is at least twice as large as KM value
  % % (consider optimal values for this criterion)
  % % log metabolite concentrations, minimum across all samples
  % c_min_over_samples = min(exp(optimal.X),[],2);
  % C_min_mat = repmat(c_min_over_samples',nr,1);
  % mark_indices = find(C_min_mat(q_info.KM_matrix_indices) > 2 * exp(qall_optimal(index.KM)));

  res_true.KM = scatter_plot(qall_true(index.KM),qall_optimal(index.KM),qall_names(index.KM),km_color,'K_{M} values [mM] (true)', 'K_{M} values [mM] (fit)', [], show_graphics,14,filenames.KM);

  switch cmb_options.parameterisation,
    case 'Keq_KV_KM_KA_KI',
      
      % figure(15); 
      res_true.Keq = scatter_plot(qall_true(index.Keq),qall_optimal(index.Keq),qall_names(index.Keq),keq_color,'K_{eq} values [unitless] (true)', 'K_{eq} values [unitless] (fit)', [], show_graphics,15,filenames.Keq);

      % figure(16); 
      res_true.Kcatf = scatter_plot(qall_true(index.Kcatf),qall_optimal(index.Kcatf),qall_names(index.Kcatf),kcat_color,'k_{cat}^{+} values [1/s] (true)', 'k_{cat}^{+} values [1/s] (fit)', [], show_graphics,16,filenames.Kcatf);
      
      % figure(17); 
      res_true.Kcatr = scatter_plot(qall_true(index.Kcatr),qall_optimal(index.Kcatr),qall_names(index.Kcatr),kcat_color,'k_{cat}^{-} values [1/s] (true)', 'k_{cat}^{-} values [1/s] (fit)', [], show_graphics,17,filenames.Kcatr);

      % figure(18); 
      res_true.Kappf = scatter_plot(qall_true(index.Kcatf),log(kapp_max.forward),qall_names(index.Kcatf),kapp_color,'k_{cat}^{+} values [1/s] (true)', 'k_{app,max}^{+} values [1/s] (fit)', [], show_graphics,18,filenames.Kappmaxf);

      % figure(19); 
      %% only if some reverse kapp values have been determined
      if sum(isfinite(kapp_max.reverse)),
        res_true.Kappr = scatter_plot(qall_true(index.Kcatr),log(kapp_max.reverse),qall_names(index.Kcatr),kapp_color,'k_{cat}^{-} values [1/s] (true)', 'k_{app,max}^{-} values [1/s] (fit)', [], show_graphics,19,filenames.Kappmaxr);
      elseif show_graphics,
        figure(19); 
        title('No reverse Kapp values determined');
        print(filenames.Kappmaxr, '-f19', '-depsc'); 
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
    case 'only_Keq'
      kin_color = pred_color;
      keq_color = fit_color;
    case 'none'
      kin_color = pred_color;
      keq_color = pred_color;
  end
  
  qall_true = cmb_q_to_qall(true.q,q_info);

  % figure(21); 
  res_true_data.c = scatter_plot(true.X,data.X.mean,[],met_color,'Metabolite levels [mM] (true)','Metabolite levels [mM] (data)', [], show_graphics,21,filenames.x);
  
  % figure(22); 
  res_true_data.e = scatter_plot(log(true.E),log(data.E.mean),[],enz_color,'Enzyme levels [mM] (true)','Enzyme levels [mM] (data)', [], show_graphics,22,filenames.e);

  % figure(23); 
  res_true_data.KV = scatter_plot(qall_true(index.KV),qall_data(index.KV),qall_names(index.KV),kin_color,'K_{V} values [1/s] (true)', 'K_{V} values [1/s] (data)', [], show_graphics,23,filenames.KV);
  
  % figure(24); 
  % % mark KM value if the minimum concentration, over all samples, is at least twice as large as KM value
  % % (consider true values for this criterion)
  % % log metabolite concentrations, minimum across all samples
  % c_min_over_samples = min(exp(true.X),[],2);
  % C_min_mat = repmat(c_min_over_samples',nr,1);
  % mark_indices = find(C_min_mat(q_info.KM_matrix_indices) > 2 * exp(qall_true(index.KM)));

  res_true_data.KM = scatter_plot(qall_true(index.KM),qall_data(index.KM),qall_names(index.KM),kin_color,'K_{M} values [mM] (true)', 'K_{M} values [mM] (data)', [], show_graphics,24,filenames.KM);

  switch cmb_options.parameterisation,
    case 'Keq_KV_KM_KA_KI',
      
      % figure(25); 
      res_true_data.Keq = scatter_plot(qall_true(index.Keq),qall_data(index.Keq),qall_names(index.Keq),keq_color,'K_{eq} values [unitless] (true)', 'K_{eq} values [unitless] (data)', [], show_graphics,25,filenames.Keq);

      % figure(26); 
      res_true_data.Kcatf = scatter_plot(qall_true(index.Kcatf),qall_data(index.Kcatf),qall_names(index.Kcatf),kin_color,'k_{cat}^{+} values [1/s] (true)', 'k_{cat}^{+} values [1/s] (data)', [], show_graphics,26,filenames.Kcatf);
      
      % figure(27); 
      res_true_data.Kcatr = scatter_plot(qall_true(index.Kcatr),qall_data(index.Kcatr),qall_names(index.Kcatr),kin_color,'k_{cat}^{-} values [1/s] (true)', 'k_{cat}^{-} values [1/s] (data)', [], show_graphics,27,filenames.Kcatr);
  
  end

end
end

% ================================================

function out = scatter_plot(xvalues,yvalues,names,color,xtitle,ytitle,mark_indices, show_graphics, fignum, filename)

eval(default('mark_indices','[]'));

ind_ok  = find(isfinite(xvalues).*isfinite(yvalues));
xvalues = xvalues(ind_ok); 
yvalues = yvalues(ind_ok); 

% all data are assumed to be on natural log scale  

geom_deviation = exp(sqrt(nanmean([yvalues(:) - xvalues(:)].^2)));

ind = find(isfinite(yvalues) .* isfinite(xvalues));
if length(ind)>2,
  linear_correlation = corrcoef(yvalues(ind),xvalues(ind)); 
  linear_correlation = linear_correlation(1,2);
else 
  linear_correlation = nan;
end

mmin = min([exp(xvalues); exp(yvalues)]);
mmax = max([exp(xvalues); exp(yvalues)]);

mmin = 10.^floor(log10(mmin));
mmax = 10.^ceil(log10(mmax));

if show_graphics,
  figure(fignum); 
  clf;
  plot([mmin,mmax],[mmin,mmax],'-','Color',[.5 .5 .5]); hold on
  plot(exp(xvalues),exp(yvalues),'.','Markersize',20,'Color',color);
  plot(exp(xvalues(mark_indices)),exp(yvalues(mark_indices)),'o','Markersize',10,'Color',color);
  %if length(names), text(exp(xvalues),exp(yvalues),names); end
  
  axis equal; axis square; axis tight
  set(gca,'XScale','log','YScale','log'); set(gca,'fontsize',16);
  title(sprintf('GeomDev: %2.2f  CorrCoeff: %2.3f',geom_deviation,linear_correlation),'FontSize',10);
  xlabel(xtitle); ylabel(ytitle);
  
  a = axis; 
  XTick = 10.^[floor(log10(min(axis))):ceil(log10(max(axis)))];
  XTickLabels = cellstr(num2str(round(log10(XTick(:))), '10^{%d}'));
  set(gca,'xtick',XTick,'ytick',XTick,'xticklabels',XTickLabels,'yticklabels',XTickLabels);

  print(filename, ['-f' num2str(fignum)], '-depsc'); 
end

out.geometric_deviation = geom_deviation;
out.linear_correlation  = linear_correlation;
