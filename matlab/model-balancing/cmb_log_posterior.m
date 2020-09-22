function [log_posterior,log_posterior_gradient] = cmb_log_posterior(y,pp,preposterior,V,cmb_options,q_info,verbose)

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

% extract X matrix and log kinetic constants from optimization vector z
[q,X] = cmb_y_to_qX(y,nm,ns);

pp.network.kinetics = cmb_q_to_kinetics(q,pp.network,cmb_options,q_info);

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

E(find(Aforward<0)) = inf;

q_log_preposterior = - 0.5 *          [q - preposterior.q.mean]' * preposterior.q.cov_inv * [q - preposterior.q.mean];
x_log_posterior    = - 0.5 * sum(sum([[X - preposterior.X.mean] ./ preposterior.X.std ].^2));
e_log_posterior    = - 0.5 * sum(sum([[E - preposterior.E.mean] ./ preposterior.E.std ].^2));

log_posterior = q_log_preposterior + x_log_posterior + e_log_posterior;

if verbose,
  q_log_preposterior
  x_log_posterior
  e_log_posterior
end


% ---------------------------------------------------------------------------
% Compute gradient

if cmb_options.use_gradient,

  if ~strcmp(cmb_options.parameterisation, 'Keq_KV_KM_KA_KI'),
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
    e_mean = preposterior.E.mean(:,it);
    e_std  = preposterior.E.std(:,it);
    log_preposterior_gradient_eX(:,it) = [[e - e_mean] ./ [e_std.^2]]' * de_dx;
    log_preposterior_gradient_eq(:,it) = [[e - e_mean] ./ [e_std.^2]]' * de_dq;
  end

  log_posterior_gradient_q = - preposterior.q.cov_inv   * [q-preposterior.q.mean] + sum(log_preposterior_gradient_eq')';
  log_posterior_gradient_X = - 1./[preposterior.X.std.^2] .* [X-preposterior.X.mean] + log_preposterior_gradient_eX;
  log_posterior_gradient   = cmb_qX_to_y(log_posterior_gradient_q,log_posterior_gradient_X,nm,ns);
  
end
