function [score,gradient] = cmb_objective(y,pp,preposterior,V,cmb_options,q_info,verbose)

% [score, gradient] = cmb_objective(y, pp, preposterior, V, cmb_options, q_info, verbose)
%
% Objective function for the optimisation with fmincon (called in cmb_estimation)
% Score function (option cmb_options.score):
%   'neg_log_posterior'      negative log posterior
%   'log_neg_log_posterior'  logarithm of negative log posterior

switch cmb_options.score,
  case 'neg_log_posterior',
    switch nargout,
      case {0,1},
        log_posterior = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,verbose);
        score         = -log_posterior; 
      case 2,
        [log_posterior,log_posterior_gradient] = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,verbose);
        score    =  log(-log_posterior); 
        gradient = -log_posterior_gradient;
    end

  case 'log_neg_log_posterior',
    switch nargout,
      case {0,1},
        log_posterior = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,verbose);
        score = log(10^-20 - log_posterior); 
      case 2,
        [log_posterior,log_posterior_gradient] = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,verbose);
        score    = log(10^-20 - log_posterior); 
        gradient = [-log_posterior_gradient] / [-log_posterior + 10^-20];
    end

  otherwise
    error('unknown score function')
end
