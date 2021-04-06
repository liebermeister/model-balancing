function preposterior = cmb_prepare_posterior(prior, data, cmb_options, q_info);

% --------------------------------------------------------
% Compute posterior distributions for X and q (formulae, see Liebermeister 2006 "Bringing .. data integration")

preposterior.X.mean = [];
preposterior.X.std  = 1 ./ sqrt( 1./[data.X.std .^2] + 1./[prior.X.std .^2] );
preposterior.X.mean = preposterior.X.std.^2 .* [ data.X.mean ./ data.X.std.^2 + prior.X.mean ./ prior.X.std.^2 ];

preposterior.V.mean = [];
preposterior.V.std  = 1 ./ sqrt( 1./[data.V.std .^2] + 1./[prior.V.std .^2] );
preposterior.V.mean = preposterior.V.std.^2 .* [ data.V.mean ./ data.V.std.^2 + prior.V.mean ./ prior.V.std.^2 ];

%preposterior.E.mean = [];
%preposterior.E.std  = 1 ./ sqrt( 1./[data.E.std .^2] + 1./[prior.E.std .^2] );
%preposterior.E.mean = preposterior.E.std.^2 .* [ data.E.mean ./ data.E.std.^2 + prior.E.mean ./ prior.E.std.^2 ];

preposterior.lnE.mean = [];
preposterior.lnE.std  = 1 ./ sqrt( 1./[data.lnE.std .^2] + 1./[prior.lnE.std .^2] );
preposterior.lnE.mean = preposterior.lnE.std.^2 .* [ data.lnE.mean ./ data.lnE.std.^2 + prior.lnE.mean ./ prior.lnE.std.^2 ];

% posterior, where data values are undefined
preposterior.X.mean(isnan(preposterior.X.mean)) = prior.X.mean(isnan(preposterior.X.mean));
preposterior.X.std(isnan(preposterior.X.std))   = prior.X.std(isnan(preposterior.X.std));

preposterior.V.mean(isnan(preposterior.V.mean)) = prior.V.mean(isnan(preposterior.V.mean));
preposterior.V.std(isnan(preposterior.V.std))   = prior.V.std(isnan(preposterior.V.std));

%preposterior.E.mean(isnan(preposterior.E.mean)) = prior.E.mean(isnan(preposterior.E.mean));
%preposterior.E.std(isnan(preposterior.E.std))   = prior.E.std(isnan(preposterior.E.std));

preposterior.lnE.mean(isnan(preposterior.lnE.mean)) = prior.lnE.mean(isnan(preposterior.lnE.mean));
preposterior.lnE.std(isnan(preposterior.lnE.std))   = prior.lnE.std(isnan(preposterior.lnE.std));

%% infer multivariate gaussian for q from multivariate gaussian for qall 
%% (e.g., from Eq 10 in Liebermeister 2006, "Bringing .. : Integration ...")

% replace nan values in data.qall.std by finite values for later calculations
data.qall.std(isnan(data.qall.std)) = nanmean(data.qall.std);

switch cmb_options.use_kinetic_data,
  
  case 'all',
    display('Using kinetic constant data');

    qall_mean           = data.qall.mean;
    qall_cov_inv        = diag(1./[data.qall.std.^2]);
    R                   = q_info.M_q_to_qall;
    %% use only non-nan data values
    ind_ok                 = find(isfinite(qall_mean));
    qall_mean              = qall_mean(ind_ok);
    qall_cov_inv           = qall_cov_inv(ind_ok,ind_ok);
    R                      = R(ind_ok,:);
    %% compute preposterior mean and covariance
    q_cov_inv              = prior.q.cov_inv + R' * qall_cov_inv * R;
    q_mean                 = q_cov_inv \ [R' * qall_cov_inv * qall_mean + prior.q.cov_inv * prior.q.mean];
    preposterior.q.mean    = q_mean;
    preposterior.q.cov_inv = q_cov_inv;
  
  case 'only_Keq',
    display('Using only equilibrium constant data');

    qall_mean           = data.qall.mean(q_info.qall.index.Keq);
    %% use narrow standard deviation for the given equilibrium constants
    qall_cov_inv        = diag(1./[0.1*ones(size(qall_mean))].^2);
    %qall_cov_inv       = diag(1./[data.qall.std(q_info.qall.index.Keq).^2]);
    R                   = q_info.M_q_to_qall(q_info.qall.index.Keq,:);
    %% use only non-nan data values
    ind_ok                 = find(isfinite(qall_mean));
    qall_mean              = qall_mean(ind_ok);
    qall_cov_inv           = qall_cov_inv(ind_ok,ind_ok);
    R                      = R(ind_ok,:);
    %% compute preposterior mean and covariance
    q_cov_inv              = prior.q.cov_inv + R' * qall_cov_inv * R;
    q_mean                 = q_cov_inv \ [R' * qall_cov_inv * qall_mean + prior.q.cov_inv * prior.q.mean];
    preposterior.q.mean    = q_mean;
    preposterior.q.cov_inv = q_cov_inv;
  
  case 'none',
    display('Using no kinetic data');
    preposterior.q = prior.q;
    
  otherwise,
    error('incorrect option');

end
