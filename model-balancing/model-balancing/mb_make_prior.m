function prior = mb_make_prior(network,q_info,mb_options);

% prior = mb_make_prior(network, q_info, mb_options)
%
% Construct prior distributions for model variables in the model balancing problem
% (from the information in mb_options.quantities, see mb_default_options)
%
% Output variables
%   prior  struct with fields
%     .q            prior for independent log kinetic constants (fields .names, .mean, .std, .prec)
%     .qdep_pseudo  pseudo prior for dependent log kinetic constants (fields .names, .mean, .std, .prec)
%     .X            prior for log metabolite concentrations (fields .mean, .std)
%     .lnE          prior for log enzyme concentrations (fields .mean, .std)
%     .V            prior for fluxes (fields .mean, .std)

% the following distribution can be improved (using the results of kinetic parameter balancing) 
% the variables for KM, KA, KI are vector containing the relevant components

[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);


% --------------------------------------------------------------
% Priors for model parameters

prior.q.names = q_info.q.names;

switch mb_options.parameterisation,
  case 'KV_KM_KA_KI',
    prior.q.mean = [mb_options.quantities.KV.mean_ln * repmat(1,nr,1); ...
                    mb_options.quantities.KM.mean_ln * repmat(1,nKM,1);  ...
                    mb_options.quantities.KA.mean_ln * repmat(1,nKA,1); ...
                    mb_options.quantities.KI.mean_ln * repmat(1,nKI,1)];
    
    prior.q.std  = [mb_options.quantities.KV.std_ln  * repmat(1,nr,1); ...
                    mb_options.quantities.KM.std_ln  * repmat(1,nKM,1); ...
                    mb_options.quantities.KA.std_ln  * repmat(1,nKA,1); ...
                    mb_options.quantities.KI.std_ln  * repmat(1,nKI,1)];
    
  case 'Keq_KV_KM_KA_KI',
    nKeqind = length(q_info.q.index.Keq);
    prior.q.mean = [mb_options.quantities.Keq.mean_ln * repmat(1,nKeqind,1); ...
                    mb_options.quantities.KV.mean_ln  * repmat(1,nr,1); ...
                    mb_options.quantities.KM.mean_ln  * repmat(1,nKM,1);  ...
                    mb_options.quantities.KA.mean_ln  * repmat(1,nKA,1); ...
                    mb_options.quantities.KI.mean_ln  * repmat(1,nKI,1)];
    
    prior.q.std  = [mb_options.quantities.Keq.std_ln * repmat(1,nKeqind,1); ...
                    mb_options.quantities.KV.std_ln  * repmat(1,nr,1); ...
                    mb_options.quantities.KM.std_ln  * repmat(1,nKM,1); ...
                    mb_options.quantities.KA.std_ln  * repmat(1,nKA,1); ...
                    mb_options.quantities.KI.std_ln  * repmat(1,nKI,1)];

  otherwise error('mb_options.parameterisation not supported');
end

prior.q.prec = diag(1./prior.q.std.^2);

% --------------------------------------------------------------
% Pseudo values for dependent model parameters

prior.qdep_pseudo.names = q_info.qdep.names;

prior.qdep_pseudo.mean = [...
    mb_options.quantities.Keq.mean_ln   * repmat(1,nr,1); ...
    mb_options.quantities.Kcatf.mean_ln * repmat(1,nr,1);  ...
    mb_options.quantities.Kcatr.mean_ln * repmat(1,nr,1)];

prior.qdep_pseudo.std  = [...
    mb_options.quantities.Keq.std_ln   * repmat(1,nr,1); ...
    mb_options.quantities.Kcatf.std_ln * repmat(1,nr,1); ...
    mb_options.quantities.Kcatr.std_ln * repmat(1,nr,1)];

prior.qdep_pseudo.prec = diag(1./prior.qdep_pseudo.std.^2);


% --------------------------------------------------------------
% Priors for metabolic variables

ns = mb_options.ns; 

prior.X.mean = log(mb_options.metabolic_prior_c_geom_mean * ones(nm, ns));
prior.lnE.mean = log(mb_options.metabolic_prior_e_geom_mean * ones(nr, ns));
prior.V.mean = zeros(nr, ns);

prior.X.std  = log(mb_options.metabolic_prior_c_geom_std) * ones(size(prior.X.mean));
prior.lnE.std = log(mb_options.metabolic_prior_e_geom_std) * ones(size(prior.lnE.mean));
prior.V.std  = 1 * ones(size(prior.V.mean)); % not used
