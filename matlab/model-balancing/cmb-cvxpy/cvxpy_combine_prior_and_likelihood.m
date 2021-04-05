function [posterior_mean, posterior_cov, posterior_std] = combine_prior_and_likelihood(prior_mean,prior_cov,data_mean,data_std)
 
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
is_ok = find(isfinite(data_mean).*isfinite(data_std));
I = eye(length(data_mean));
P = I(is_ok,:);

data_mean = data_mean(is_ok);
data_std  = data_std(is_ok);
  
prior_cov_inv     = inv(prior_cov);
data_cov_inv      = diag(1./data_std.^2);
posterior_cov_inv = prior_cov_inv + P' * data_cov_inv * P;
posterior_cov     = inv(posterior_cov_inv);
posterior_std     = sqrt(diag(posterior_cov));
posterior_mean    = posterior_cov_inv \ [prior_cov_inv * prior_mean + P' * data_cov_inv * data_mean];