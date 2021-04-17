function [log_posterior,log_posterior_gradient] = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,verbose)

% [log_posterior,log_posterior_gradient] = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,verbose)

% ------------------------------
% Initialise some variables
  
eval(default('verbose','0'))  

no_warnings = 1;

global LP_info % variable defined in cmb_estimation

eval(default('verbose','0'));

if verbose,
  if sum(LP_info.y_ineq_A * y > LP_info.y_ineq_b - LP_info.epsilon) ~=0,
    warning('Constraint violation detected');
    max_violation = max(LP_info.y_ineq_A * y - [ LP_info.y_ineq_b - LP_info.epsilon])
  end
end
  
[nr,nm,nx,KM_indices,KA_indices,KI_indices,nKM,nKA,nKI] = network_numbers(pp.network);

ns = size(preposterior.X.mean,2);


% ------------------------------
% extract X matrix and log kinetic constants from optimization vector z

[q,X] = cmb_y_to_qX(y,nm,ns);

pp.network.kinetics = cmb_q_to_kinetics(q,pp.network,cmb_options,q_info);

q_log_preposterior = - 0.5 *          [q - preposterior.q.mean]' * preposterior.q.cov_inv * [q - preposterior.q.mean];

x_log_posterior    = - 0.5 * sum(sum([[X - preposterior.X.mean] ./ preposterior.X.std ].^2));


% ------------------------------
% compute enzyme levels (ECM code)

for it = 1:ns,
  v      = V(:,it);
  x      = X(:,it);
  c      = exp(x);
  pp.v   = v;
  [~, e] = ecm_get_score(cmb_options.ecm_score,x,pp,no_warnings);
  E(:,it) = e;
  Aforward(:,it) = RT * diag(sign(v)) * [log(pp.network.kinetics.Keq) - pp.network.N' * x];
end

if verbose,
  if find(Aforward<0), 
    display('Constraint violation detected');
  end
end

% Avoid zero enzyme levels (because of description on log scale)
E(find(V==0)) = 10^-10;
E(find([Aforward<0])) = inf;

ln_e_log_posterior_upper  = [[log(E) - preposterior.lnE.mean] ./ preposterior.lnE.std ].^2 .* double(log(E) >= preposterior.lnE.mean);
ln_e_log_posterior_lower  = [[log(E) - preposterior.lnE.mean] ./ preposterior.lnE.std ].^2 .* double(log(E) < preposterior.lnE.mean);

switch cmb_options.enzyme_score_type,
  case 'quadratic'
    % normal quadratic term: with this option, MB is not guaranteed to be convex!
    ln_e_log_posterior  = - 0.5 * sum(sum(ln_e_log_posterior_lower + ln_e_log_posterior_upper));
  case 'monotonic',
    % set likehood for E values below posterior mean to 0 
    % -> makes enzyme posterior score monotonically increasing;
    % and therefore overall MP problem convex
    ln_e_log_posterior = - 0.5 * sum(sum(ln_e_log_posterior_upper));
  case 'interpolated',
    ln_e_log_posterior = - 0.5 * sum(sum(cmb_options.enzyme_score_alpha * ln_e_log_posterior_lower + ln_e_log_posterior_upper));
    otherwise('error');
end

% -----------------------------------------------
% c / KM pseudo values term

ln_c_over_km_log_posterior = 0;

if cmb_options.beta_ln_c_over_km > 0,
  for it = 1:size(X,2),  
    x = X(:,it);
    dum = [repmat(x',nr,1) - log(pp.network.kinetics.KM)];
    dum = sum(find(pp.network.kinetics.KM~=0));
    c_over_km_score = 1/2 * cmb_options.beta_ln_c_over_km^2 * sum(sum(dum.^2));
    ln_c_over_km_log_posterior = ln_c_over_km_log_posterior + c_over_km_score;
  end
% sigma for pseudo value term for sum_til ln(c_i(t)/km_il)
end

log_posterior = q_log_preposterior + x_log_posterior + ln_e_log_posterior + ln_c_over_km_log_posterior;

if verbose,
  q_log_preposterior
  x_log_posterior
  ln_e_log_posterior
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
    %e_mean = preposterior.E.mean(:,it);
    %e_std  = preposterior.E.std(:,it);
    ln_e_mean = preposterior.lnE.mean(:,it);
    ln_e_std  = preposterior.lnE.std(:,it);
    %log_preposterior_gradient_eX(:,it) = [[e - e_mean] ./ [e_std.^2]]' * de_dx;
    %log_preposterior_gradient_eq(:,it) = [[e - e_mean] ./ [e_std.^2]]' * de_dq;
    log_preposterior_gradient_eX(:,it) = [[log(e) - ln_e_mean] ./ [ln_e_std.^2]]' * de_dx;
    log_preposterior_gradient_eq(:,it) = [[log(e) - ln_e_mean] ./ [ln_e_std.^2]]' * de_dq;
  end

  log_posterior_gradient_q = - preposterior.q.cov_inv * [q-preposterior.q.mean] + sum(log_preposterior_gradient_eq')';
  log_posterior_gradient_X = - [X-preposterior.X.mean]./[preposterior.X.std.^2] + log_preposterior_gradient_eX;
  log_posterior_gradient   = cmb_qX_to_y(log_posterior_gradient_q,log_posterior_gradient_X,nm,ns);
end
