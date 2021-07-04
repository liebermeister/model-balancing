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

index       = q_info.qall.index;
M_q_to_qall = q_info.M_q_to_qall;
M_q_to_qdep = q_info.M_q_to_qdep;

switch cmb_options.use_kinetic_data,
  case 'all',
    display('Using kinetic and equilibrium constants data');
    qall_data_mean = data.qall.mean;
    qall_data_std  = data.qall.std;    
    M_q_to_qdata   = q_info.M_q_to_qall;

  case 'only_Keq',
    display('Using equilibrium constants data');
    qall_data_mean = data.qall.mean(index.Keq);
    qall_data_std  = data.qall.std(index.Keq);
    M_q_to_qdata   = q_info.M_q_to_qall(index.Keq,:);
  
  case 'none',
    display('cmb_prepare_posterior: Not using kinetic or equilibrium constants data');
    qall_data_mean = [];
    qall_data_std  = [];
    M_q_to_qdata   = [];
    
  otherwise,
    error('incorrect option');

end

%% use only non-nan data values
I               = eye(length(data.qall.mean));
is_ok           = find(isfinite(qall_data_mean) .* isfinite(qall_data_std));
M_qall_to_qdata = I(is_ok,:);
M_q_to_qdata    = M_q_to_qdata(is_ok,:);
qdata_mean      = qall_data_mean(is_ok);
qdata_std       = qall_data_std(is_ok);
qdata_prec      = diag(1./[qdata_std.^2]);

% compute preposterior mean and covariance

% OLD VARIANT: compute correlated posterior for q (from prior and data)

if length(qdata_mean),
  if cmb_options.use_pseudo_values,
    %%  with pseudo values
    q_prec = prior.q.prec + M_q_to_qdata' * qdata_prec * M_q_to_qdata + M_q_to_qdep' * prior.qdep_pseudo.prec * M_q_to_qdep;
    q_mean = q_prec \ [M_q_to_qdata' * qdata_prec * qdata_mean + prior.q.prec * prior.q.mean + M_q_to_qdep' * prior.qdep_pseudo.prec * prior.qdep_pseudo.mean];
  else,
    %% without pseudo values
    q_prec = prior.q.prec + M_q_to_qdata' * qdata_prec * M_q_to_qdata;
    q_mean = q_prec \ [M_q_to_qdata' * qdata_prec * qdata_mean + prior.q.prec * prior.q.mean];
  end
else
  q_prec = prior.q.prec;
  q_mean = prior.q.mean;
end

% preposterior for basic kinetic constants

preposterior.q.names   = q_info.q.names;
preposterior.q.mean    = q_mean;
preposterior.q.prec    = q_prec;
preposterior.q.cov     = inv(q_prec);

% preposterior for all kinetic constants

preposterior.qall.names= q_info.qall.names;
preposterior.qall.mean = M_q_to_qall * q_mean;
preposterior.qall.cov  = M_q_to_qall * preposterior.q.cov * M_q_to_qall'; 
preposterior.qall.prec = pinv(preposterior.qall.cov);

% preposterior for individual types of kinetic constants (for variant "use_crosscovariances =0")

preposterior.qtypes.Keq.mean   = preposterior.qall.mean(index.Keq);
preposterior.qtypes.KM.mean    = preposterior.qall.mean(index.KM);
preposterior.qtypes.KA.mean    = preposterior.qall.mean(index.KA);
preposterior.qtypes.KI.mean    = preposterior.qall.mean(index.KI); 
preposterior.qtypes.Kcatf.mean = preposterior.qall.mean(index.Kcatf);
preposterior.qtypes.Kcatr.mean = preposterior.qall.mean(index.Kcatr);
preposterior.qtypes.Keq.prec   = inv(preposterior.qall.cov(index.Keq,index.Keq));
preposterior.qtypes.KM.prec    = inv(preposterior.qall.cov(index.KM,index.KM));
preposterior.qtypes.KA.prec    = inv(preposterior.qall.cov(index.KA,index.KA));
preposterior.qtypes.KI.prec    = inv(preposterior.qall.cov(index.KI,index.KI)); 
preposterior.qtypes.Kcatf.prec = inv(preposterior.qall.cov(index.Kcatf,index.Kcatf));
preposterior.qtypes.Kcatr.prec = inv(preposterior.qall.cov(index.Kcatr,index.Kcatr));



% % NEW variant: compute everything for qall, then select vectors and matrices for q
% 
% % combine prior + pseudo values -> pseudoprior for all constants
% M_q_to_qall      = q_info.M_q_to_qall;
% M_qdep_to_qall   = q_info.M_qdep_to_qall;
% qall_pprior_prec = M_q_to_qall * prior.q.prec * M_q_to_qall' + M_qdep_to_qall * prior.qdep_pseudo.prec * M_qdep_to_qall';
% qall_pprior_mean = qall_pprior_prec \ [M_q_to_qall * prior.q.prec * prior.q.mean + M_qdep_to_qall * prior.qdep_pseudo.prec * prior.qdep_pseudo.mean];
% 
% % preposterior mean and covariance
% if length(qdata_mean),
%   qall_prec = qall_pprior_prec + M_qall_to_qdata' * qdata_prec * M_qall_to_qdata;
%   qall_mean = qall_prec \ [M_qall_to_qdata' * qdata_prec * qdata_mean + qall_pprior_prec * qall_pprior_mean];
% else
%   % important to write this explicitly, matlab has a strange bug when qdata_mean is an empty vector
%   qall_prec = qall_pprior_prec;
%   qall_mean = qall_pprior_mean;
% end
% 
% preposterior.qall.mean = qall_mean;
% preposterior.qall.prec = qall_prec;
% preposterior.qall.cov  = inv(qall_prec);
% 
% preposterior.q.mean    = qall_mean(q_info.q.index_q);
% preposterior.q.cov     = preposterior.qall.cov(q_info.q.index_q,q_info.q.index_q);
% preposterior.q.prec    = pinv(preposterior.q.cov);
