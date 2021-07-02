function [log_posterior,log_posterior_gradient] = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,verbose)

% [log_posterior,log_posterior_gradient] = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,verbose)

% ------------------------------
% Initialise some variables
  
eval(default('verbose','0'))  

global LP_info % variable defined in cmb_estimation

if verbose,
  if length(LP_info),
    if sum(LP_info.y_ineq_A * y > LP_info.y_ineq_b - LP_info.epsilon) ~=0,
      warning('Constraint violation detected');
      max_violation = max(LP_info.y_ineq_A * y - [ LP_info.y_ineq_b - LP_info.epsilon])
    end
  end
end

no_warnings = 1;

[nr,nm,nx,KM_indices,KA_indices,KI_indices,nKM,nKA,nKI] = network_numbers(pp.network);

ns = size(preposterior.X.mean,2);

% ------------------------------
% Extract X matrix and log kinetic constants from optimization vector z

[q,X] = cmb_y_to_qX(y,nm,ns);

% version A: consider only basic kinetic constants  
log_posterior_q = - 0.5 * [q - preposterior.q.mean]' * preposterior.q.prec * [q - preposterior.q.mean];

% äquivalent zu version A
% consider only basic kinetic constants
% mean and prec for "qmod" (vector of kinetic constants as they are used in the model)
% index           = q_info.qall.index;
% ind             = [index.Keq, index.KM, index.KA, index.KI, index.Kcatf, index.Kcatr];
% M_q_qmod        = q_info.M_q_to_qall(ind,:);
% iM_q_qmod       = pinv(M_q_qmod);
% qall            = cmb_q_to_qall(q, q_info);
% qmod            = qall(ind);
% qmod_mean       = M_q_qmod * preposterior.q.mean;
% qmod_prec       = iM_q_qmod' * preposterior.q.prec * iM_q_qmod;
% log_posterior_q = - 0.5 * [qmod - qmod_mean]' * qmod_prec * [qmod - qmod_mean];

% äquivalent zu version A
% consider all kinetic constants
% index     = q_info.qall.index;
% ind       = [index.Keq, index.KM, index.KA, index.KI, index.Kcatf, index.Kcatr];
% qall      = cmb_q_to_qall(q, q_info);
% qall_mean = preposterior.qall.mean;
% qall_prec = preposterior.qall.prec;
% qall_cov  = preposterior.qall.cov; 
% log_posterior_q = - 0.5 * [qall(ind) - qall_mean(ind)]' * pinv(qall_cov(ind,ind)) * [qall(ind) - qall_mean(ind)];

% version B: consider all kinetic constants (ignore convariances between them)
% nicht äquivalent zu version A
% index     = q_info.qall.index;
% ind       = [index.Keq, index.KM, index.KA, index.KI, index.Kcatf, index.Kcatr];
% qall      = cmb_q_to_qall(q, q_info);
% qall_mean = preposterior.qall.mean;
% qall_prec = preposterior.qall.prec;
% qall_cov  = preposterior.qall.cov; 
% log_posterior_q = - 0.5 * [qall(ind) - qall_mean(ind)]' * qall_prec(ind,ind) * [qall(ind) - qall_mean(ind)];

% äquivalent zu version B
% per variable type 
% index     = q_info.qall.index;
% ind       = [index.Keq, index.KM, index.KA, index.KI, index.Kcatf, index.Kcatr];
% qall      = cmb_q_to_qall(q, q_info);
% qall_cov  = preposterior.qall.cov;
% qall_mean = preposterior.qall.mean;
% log_posterior_q = - 0.5 * [...
%       [qall(index.Keq)   - qall_mean(index.Keq)]'   * inv(qall_cov(index.Keq,index.Keq))     * [qall(index.Keq)   - qall_mean(index.Keq)] ...
%     + [qall(index.KM)    - qall_mean(index.KM)]'    * inv(qall_cov(index.KM,index.KM))       * [qall(index.KM)    - qall_mean(index.KM)] ...
%     + [qall(index.KA)    - qall_mean(index.KA)]'    * inv(qall_cov(index.KA,index.KA))       * [qall(index.KA)    - qall_mean(index.KA)] ...
%     + [qall(index.KI)    - qall_mean(index.KI)]'    * inv(qall_cov(index.KI,index.KI))       * [qall(index.KI)    - qall_mean(index.KI)] ...
%     + [qall(index.Kcatf) - qall_mean(index.Kcatf)]' * inv(qall_cov(index.Kcatf,index.Kcatf)) * [qall(index.Kcatf) - qall_mean(index.Kcatf)] ...
%     + [qall(index.Kcatr) - qall_mean(index.Kcatr)]' * inv(qall_cov(index.Kcatr,index.Kcatr)) * [qall(index.Kcatr) - qall_mean(index.Kcatr)]];

% äquivalent zu version B
% per variable type 
% index     = q_info.qall.index;
% ind       = [index.Keq, index.KM, index.KA, index.KI, index.Kcatf, index.Kcatr];
% qall      = cmb_q_to_qall(q, q_info);
% qall_cov  = preposterior.qall.cov;
% qall_prec = preposterior.qall.prec;
% qall_mean = preposterior.qall.mean;
% log_posterior_q = - 0.5 * [...
%       [qall(index.Keq)   - qall_mean(index.Keq)]'   * qall_prec(index.Keq,index.Keq)     * [qall(index.Keq)   - qall_mean(index.Keq)] ...
%     + [qall(index.KM)    - qall_mean(index.KM)]'    * qall_prec(index.KM,index.KM)       * [qall(index.KM)    - qall_mean(index.KM)] ...
%     + [qall(index.KA)    - qall_mean(index.KA)]'    * qall_prec(index.KA,index.KA)       * [qall(index.KA)    - qall_mean(index.KA)] ...
%     + [qall(index.KI)    - qall_mean(index.KI)]'    * qall_prec(index.KI,index.KI)       * [qall(index.KI)    - qall_mean(index.KI)] ...
%     + [qall(index.Kcatf) - qall_mean(index.Kcatf)]' * qall_prec(index.Kcatf,index.Kcatf) * [qall(index.Kcatf) - qall_mean(index.Kcatf)] ...
%     + [qall(index.Kcatr) - qall_mean(index.Kcatr)]' * qall_prec(index.Kcatr,index.Kcatr) * [qall(index.Kcatr) - qall_mean(index.Kcatr)]];

% ------------------------------

log_posterior_x = - 0.5 * sum(sum([[X - preposterior.X.mean] ./ preposterior.X.std ].^2));

% ------------------------------
% compute enzyme levels (ECM code)

pp.network.kinetics = cmb_q_to_kinetics(q,pp.network,cmb_options,q_info);

for it = 1:ns,
  v      = V(:,it);
  x      = X(:,it);
  c      = exp(x);
  pp.v   = v;
  [~, E] = ecm_get_score(cmb_options.ecm_score,x,pp,no_warnings);
  lnE(:,it) = log(E);
  A(:,it) = RT * [log(pp.network.kinetics.Keq) - pp.network.N' * x];
end

% check
if verbose,
  if find(V.*A<0),
    display('WARNING: Constraint violation detected!');
  end
end

% Avoid zero or negative enzyme levels (because of description on log scale)
lnE(find(V==0))   = log(10^-10);
lnE(find(V.*A<0)) = inf;

lnE_log_posterior_upper  = - 0.5 * sum(sum( [[lnE - preposterior.lnE.mean] ./ preposterior.lnE.std ].^2 .* double(lnE >= preposterior.lnE.mean)));
lnE_log_posterior_lower  = - 0.5 * sum(sum( [[lnE - preposterior.lnE.mean] ./ preposterior.lnE.std ].^2 .* double(lnE <  preposterior.lnE.mean)));

switch cmb_options.enzyme_score_type,
  case 'quadratic'
    % with this option, MB is not guaranteed to be convex!
    my_alpha = 1;
  case 'monotonic',
    % set likehood for E values below posterior mean to 0 
    % -> monotonically increasing enzyme posterior score; overall MP problem is convex
    my_alpha = 0;
  case 'interpolated',
    my_alpha = cmb_options.enzyme_score_alpha;
  otherwise('error');
end

log_posterior_lnE = my_alpha * lnE_log_posterior_lower + lnE_log_posterior_upper;

% -----------------------------------------------
% c / KM pseudo values term

log_posterior_ln_c_over_km = 0;

if cmb_options.beta_ln_c_over_km > 0,
  %% 1/sigma for pseudo value term for sum_til ln(c_i(t)/km_il)
  for it = 1:size(X,2),  
    x = X(:,it);
    log_c_over_km   = [repmat(x',nr,1) - log(pp.network.kinetics.KM)];
    log_c_over_km(find(pp.network.kinetics.KM==0)) = 0;
    c_over_km_score = - 0.5 * cmb_options.beta_ln_c_over_km^2 * sum(sum(log_c_over_km.^2));
    log_posterior_ln_c_over_km = log_posterior_ln_c_over_km + c_over_km_score;
  end
end

% -----------------------------------------------

log_posterior = log_posterior_q + log_posterior_x + log_posterior_lnE + log_posterior_ln_c_over_km;

if verbose,
  % show z-scores
  z_q = -2*log_posterior_q
  z_x =-2*log_posterior_x
  z_lne = -2*log_posterior_lnE
  z_lnckm = -2*log_posterior_ln_c_over_km
end


% ---------------------------------------------------------------------------
% Compute gradient

if cmb_options.use_gradient,

  switch cmb_options.enzyme_score_type,
  case 'quadratic',
  otherwise,
    error('Gradient not supported'); 
end

  if ~strcmp(cmb_options.parameterisation, 'Keq_KV_KM_KA_KI'),
    error('Gradient not supported'); 
  end

  if cmb_options.beta_ln_c_over_km ~= 0,
    error('Gradient not supported'); 
  end

  %% Compute posterior terms due to enzyme levels
  for it = 1:ns,
    v = V(:,it);
    x = X(:,it);
    e = E(:,it);
    c = exp(x);
    k = exp(q);
    %% Compute metabolite and parameter elasticties
    pp.network.kinetics.u = e;
    [Eun_c,Eun_p,parameters] = elasticities(pp.network,c);
    %% reorder parameter elasticities (different convention in elasticities.m)
    for itt = 1:nr,
      Eun_kin      = ms_vector2par(Eun_p(itt,:)',pp.network.kinetics,pp.network);
      Eun_k(itt,:) = [Eun_kin.Keq; Eun_kin.KV; Eun_kin.KM(KM_indices); Eun_kin.KA(KA_indices); Eun_kin.KI(KI_indices)]';
    end
    %% Compute posterior terms due to enzyme levels
    de_dx  = - diag(e./v) * Eun_c * diag(c);
    de_dq  = - diag(e./v) * Eun_k * diag(k);
    lnE_mean = preposterior.lnE.mean(:,it);
    lnE_std  = preposterior.lnE.std(:,it);
    log_preposterior_gradient_eX(:,it) = [[log(e) - lnE_mean] ./ [lnE_std.^2]]' * de_dx;
    log_preposterior_gradient_eq(:,it) = [[log(e) - lnE_mean] ./ [lnE_std.^2]]' * de_dq;
  end

  log_posterior_gradient_q = - preposterior.q.prec * [q-preposterior.q.mean] + sum(log_preposterior_gradient_eq')';
  log_posterior_gradient_X = - [X-preposterior.X.mean]./[preposterior.X.std.^2] + log_preposterior_gradient_eX;
  log_posterior_gradient   = cmb_qX_to_y(log_posterior_gradient_q,log_posterior_gradient_X,nm,ns);
  
end
