function [posterior_mean, posterior_cov, posterior_std, posterior_prec] = combine_prior_and_likelihood(prior_mean,prior_cov,data_mean,data_std)
 
% [posterior_mean, posterior_cov] = combine_prior_and_likelihood(prior_mean,prior_cov,data_mean,data_std)
%
% combine prior and likelihood information for a vector of variables (on log scale)
% use formula from convenience kinetics / Bayesian estimation paper
% 
% Inputs:
%   prior_mean: vector of prior means
%   prior_cov:  prior covariance matrix
%   data_mean:  data values (one value for each variable)
%   data_std:   data standard deviations (one value for each variable)
% Outputs:
%   posterior_mean: vector of posterior means
%   posterior_cov:  posterior covariance matrix

% projector matrix P (from all variables to those for which valid data values are available)

is_ok          = find(isfinite(data_mean) .* isfinite(data_std));
I              = eye(length(data_mean));
P              = I(is_ok,:);
data_mean      = data_mean(is_ok);
data_std       = data_std(is_ok);  
prior_prec     = inv(prior_cov);
data_prec      = diag(1./data_std.^2);
posterior_prec = prior_prec + P' * data_prec * P;
posterior_mean = posterior_prec \ [prior_prec * prior_mean + P' * data_prec * data_mean];
posterior_cov  = inv(posterior_prec);
posterior_std  = sqrt(diag(posterior_cov));
