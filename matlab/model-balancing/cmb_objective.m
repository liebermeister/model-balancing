function   [score,gradient] = cmb_objective(y,pp,preposterior,V,cmb_options,q_info)

%% CURRENTLY NOT USED  
  
switch cmb_options.score,
  case 'neg_log_posterior',
    switch nargout,
      case 1,
        log_posterior = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info);
        score    = -log_posterior; 
      case 2,
        [log_posterior,log_posterior_gradient] = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info);
        score    =  log(-log_posterior); 
        gradient = -log_posterior_gradient;
    end
    
  case 'log_neg_log_posterior',
    switch nargout,
      case 1,
        log_posterior         = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info);
        score = log(-log_posterior + 10^-20); 
      case 2,
        [log_posterior,log_posterior_gradient] = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info);
        score    = log(-log_posterior + 10^-20); 
        gradient = [-log_posterior_gradient] / [-log_posterior + 10^-20];
    end
end
