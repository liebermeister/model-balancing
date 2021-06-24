function preposterior = cmb_prepare_posterior(prior, data, cmb_options, q_info);

% --------------------------------------------------------
% Compute posterior distributions for X and q (formulae, see Liebermeister 2006 "Bringing .. data integration")

preposterior.X.mean = [];
preposterior.X.std  = 1 ./ sqrt( 1./[data.X.std .^2] + 1./[prior.X.std .^2] );
preposterior.X.mean = preposterior.X.std.^2 .* [ data.X.mean ./ data.X.std.^2 + prior.X.mean ./ prior.X.std.^2 ];

preposterior.V.mean = [];
preposterior.V.std  = 1 ./ sqrt( 1./[data.V.std .^2] + 1./[prior.V.std .^2] );
preposterior.V.mean = preposterior.V.std.^2 .* [ data.V.mean ./ data.V.std.^2 + prior.V.mean ./ prior.V.std.^2 ];

preposterior.lnE.mean = [];
preposterior.lnE.std  = 1 ./ sqrt( 1./[data.lnE.std .^2] + 1./[prior.lnE.std .^2] );
preposterior.lnE.mean = preposterior.lnE.std.^2 .* [ data.lnE.mean ./ data.lnE.std.^2 + prior.lnE.mean ./ prior.lnE.std.^2 ];

% posterior, where data values are undefined
preposterior.X.mean(isnan(preposterior.X.mean)) = prior.X.mean(isnan(preposterior.X.mean));
preposterior.X.std(isnan(preposterior.X.std))   = prior.X.std(isnan(preposterior.X.std));

preposterior.V.mean(isnan(preposterior.V.mean)) = prior.V.mean(isnan(preposterior.V.mean));
preposterior.V.std(isnan(preposterior.V.std))   = prior.V.std(isnan(preposterior.V.std));

preposterior.lnE.mean(isnan(preposterior.lnE.mean)) = prior.lnE.mean(isnan(preposterior.lnE.mean));
preposterior.lnE.std(isnan(preposterior.lnE.std))   = prior.lnE.std(isnan(preposterior.lnE.std));

% --------------------------------------------------------------------------------
% infer multivariate gaussian for q from multivariate gaussian for qall 
% (e.g., from Eq 10 in Liebermeister 2006, "Bringing .. : Integration ...")

% replace nan values in data.qall.std by finite values for later calculations
% data.qall.std(isnan(data.qall.std)) = nanmean(data.qall.std);

switch cmb_options.use_kinetic_data,
  case 'all',
    display('Using kinetic and equilibrium constants data');
    qall_data_mean = data.qall.mean;
    qall_data_std  = data.qall.std;    
    R              = q_info.M_q_to_qall;

  case 'only_Keq',
    display('Using equilibrium constants data');
    qall_data_mean = data.qall.mean(q_info.qall.index.Keq);
    qall_data_std  = data.qall.std(q_info.qall.index.Keq);
    R              = q_info.M_q_to_qall(q_info.qall.index.Keq,:);
  
  case 'none',
    display('cmb_prepare_posterior: Not using kinetic or equilibrium constants data');
    qall_data_mean = [];
    qall_data_std  = [];
    R              = [];
    
  otherwise,
    error('incorrect option');

end

%% use only non-nan data values
is_ok             = find(isfinite(qall_data_mean).*isfinite(qall_data_std));
I                 = eye(length(data.qall.mean));
P                 = I(is_ok,:);
qall_data_mean    = qall_data_mean(is_ok);
qall_data_std     = qall_data_std(is_ok);
qall_data_cov_inv = diag(1./[qall_data_std.^2]);
R                 = R(is_ok,:);

% compute preposterior mean and covariance

% OLD VARIANT: compute everything only for q
% q_cov_inv              = prior.q.cov_inv + R' * qall_data_cov_inv * R;
% q_mean                 = q_cov_inv \ [R' * qall_data_cov_inv * qall_data_mean + prior.q.cov_inv * prior.q.mean];
% preposterior.q.mean    = q_mean;
% preposterior.q.cov_inv = q_cov_inv;
% preposterior.q.cov     = inv(preposterior.q.cov_inv);

% NEW variant: compute everything for qall, then select vectors and matrices for q

% prior
M_q_to_qall         = q_info.M_q_to_qall;
M_qdep_to_qall      = q_info.M_qdep_to_qall;
D1 = inv(M_q_to_qall * prior.q.cov_inv * M_q_to_qall' + M_qdep_to_qall * prior.qdep_pseudo.cov_inv * M_qdep_to_qall');
D2 = D1 * [M_q_to_qall * prior.q.cov_inv * prior.q.mean + M_qdep_to_qall * prior.qdep_pseudo.cov_inv * prior.qdep_pseudo.mean];
qall_prior_mean     = D2;
qall_prior_cov      = D1;
qall_prior_std      = sqrt(diag(qall_prior_cov));
qall_prior_cov_inv  = inv(qall_prior_cov);

% preposterior mean and covariance
if length(qall_data_mean),
  qall_cov_inv        = qall_prior_cov_inv + P' * qall_data_cov_inv * P;
  qall_mean           = qall_cov_inv \ [qall_prior_cov_inv * qall_prior_mean + P' * qall_data_cov_inv * qall_data_mean];
  qall_cov            = inv(qall_cov_inv);
else
  % important to write this explicitly, matlab has a strange bug when qall_data_mean is an empty vector
  qall_mean           = qall_prior_mean;
  qall_cov_inv        = qall_prior_cov_inv;
  qall_cov            = inv(qall_prior_cov_inv);
end

preposterior.q.mean = qall_mean(q_info.q.index_q);
preposterior.q.cov  = qall_cov(q_info.q.index_q,q_info.q.index_q);
preposterior.q.cov_inv = qall_cov_inv(q_info.q.index_q,q_info.q.index_q);

%preposterior.q.mean
%preposterior.q.cov
