function gradient = mb_log_posterior_gradient(q, X, lnE, V, pp, preposterior, mb_options, q_info, alpha)

% gradient = mb_log_posterior_gradient(q, X, lnE, V, pp, preposterior, mb_options, q_info, alpha)
%
% Gradient of the log posterior score (see mb_log_posterior) with respect to the vector y
% of optimisation variables. Only for the common modular rate law (ecm_score = 'emc4cm').
%
% Enzyme levels are given by ln e = ln |v| - ln |w(x,q)|, where w is the rate per enzyme:
%   w = R * (kplus * theta_plus - kminus * theta_minus) / D
%
% The score for enzyme levels has a continuous first derivative at ln e = data mean,
% so this gradient exists for all values of alpha.
%
% Input variables
%   q, X    log kinetic constants and log metabolite concentrations
%   lnE     log enzyme levels (reactions x metabolic states), as computed in mb_log_posterior
%   V       flux matrix (reactions x metabolic states)
%   alpha   weight of the enzyme score for enzyme levels below the preposterior mean

if ~strcmp(mb_options.ecm_score,'emc4cm'),
  error('Gradient only available for ecm_score = ''emc4cm''');
end

[nm,ns] = size(X);

% ------------------------------
% Terms from the preposterior of kinetic constants and metabolite levels

if mb_options.use_crosscovariances,
  gradient_q = - preposterior.q.prec * [q - preposterior.q.mean];
else
  qall       = mb_q_to_qall(q, q_info);
  gradient_q = zeros(size(q));
  for type = {'Keq','KM','KA','KI','Kcatf','Kcatr'},
    ind        = q_info.qall.index.(type{1});
    prior_type = preposterior.qtypes.(type{1});
    gradient_q = gradient_q - q_info.M_q_to_qall(ind,:)' * prior_type.prec * [qall(ind) - prior_type.mean];
  end
end

gradient_X = - [X - preposterior.X.mean] ./ preposterior.X.std.^2;

% ------------------------------
% Term from enzyme levels: d(score)/dy = - weight * z * d(ln e)/dy = + weight * z * d(ln w)/dy

weight = ones(size(lnE));
weight(lnE < preposterior.lnE.mean) = alpha;
weighted_z = weight .* [lnE - preposterior.lnE.mean] ./ preposterior.lnE.std.^2;
weighted_z(V==0) = 0; % reactions without flux do not contribute to the score

for it = 1:ns,
  [dlnw_dx, dlnw_dq] = rate_elasticities(X(:,it), q, pp.network, q_info);
  gradient_X(:,it)   = gradient_X(:,it) + dlnw_dx' * weighted_z(:,it);
  gradient_q         = gradient_q       + dlnw_dq' * weighted_z(:,it);
end

% ------------------------------
% Term from c/KM pseudo values

beta = mb_options.beta_ln_c_over_km;

if beta > 0,
  [nr,~,~,KM_ind] = network_numbers(pp.network);
  [~,ind_m]       = ind2sub([nr,nm],KM_ind);
  ind_KM          = q_info.q.index.KM(:);
  for it = 1:ns,
    ln_c_over_km       = X(ind_m,it) - q(ind_KM);
    gradient_X(:,it)   = gradient_X(:,it) - beta^2 * accumarray(ind_m(:), ln_c_over_km, [nm,1]);
    gradient_q(ind_KM) = gradient_q(ind_KM) + beta^2 * ln_c_over_km;
  end
end

gradient = mb_qX_to_y(gradient_q, gradient_X, nm, ns);


% --------------------------------------------------------------------------------

function [dlnw_dx, dlnw_dq] = rate_elasticities(x, q, network, q_info)

% Derivatives of ln |w| (rate per enzyme, common modular rate law, see modular_velocities)
% with respect to log concentrations x (reactions x metabolites) and log constants q (reactions x q)

[nr,nm,~,KM_ind,KA_ind,KI_ind] = network_numbers(network);

h     = network.kinetics.h(:);
index = q_info.q.index;

ln_Keq = q_info.M_qKeqind_to_qallKeq * q(index.Keq);
ln_KV  = q(index.KV);
ln_KM  = q(index.KM);  ln_KM = ln_KM(:);
ln_KA  = q(index.KA);  ln_KA = ln_KA(:);
ln_KI  = q(index.KI);  ln_KI = ln_KI(:);

% Substrates and products: one element for each KM value

[r,m]      = ind2sub([nr,nm],KM_ind);
r = r(:); m = m(:);
n          = full(network.N(sub2ind([nm,nr],m,r)));
order_plus = h(r) .* max(-n,0);
order_minus= h(r) .* max( n,0);

a   = x(m) - ln_KM;               % ln(c/KM)
sat = 1 ./ (1 + exp(-a));         % c / (c + KM)

ln_theta_plus  = accumarray(r, order_plus  .* a, [nr,1]);
ln_theta_minus = accumarray(r, order_minus .* a, [nr,1]);
psi_plus       = exp(accumarray(r, order_plus  .* log1p(exp(a)), [nr,1]));
psi_minus      = exp(accumarray(r, order_minus .* log1p(exp(a)), [nr,1]));
D              = psi_plus + psi_minus - 1;

ln_haldane  = ln_Keq - accumarray(r, n .* ln_KM, [nr,1]);
term_plus   = exp(ln_KV + 0.5 * h .* ln_haldane + ln_theta_plus);
term_minus  = exp(ln_KV - 0.5 * h .* ln_haldane + ln_theta_minus);
share_plus  = term_plus  ./ [term_plus - term_minus];
share_minus = term_minus ./ [term_plus - term_minus];
share_sum   = share_plus + share_minus;

dlnw_da = share_plus(r) .* order_plus - share_minus(r) .* order_minus ...
          - [psi_plus(r) .* order_plus + psi_minus(r) .* order_minus] .* sat ./ D(r);

% Regulation: activators (factor c/(c+KA)) and inhibitors (factor KI/(c+KI))

W          = network.regulation_matrix;
[r_A,m_A]  = ind2sub([nr,nm],KA_ind);
[r_I,m_I]  = ind2sub([nr,nm],KI_ind);
r_A = r_A(:); m_A = m_A(:); r_I = r_I(:); m_I = m_I(:);
dlnw_dlnKA = - full(W(KA_ind)) ./ (1 + exp(x(m_A) - ln_KA));
dlnw_dlnKI =   abs(full(W(KI_ind))) ./ (1 + exp(ln_KI - x(m_I)));
dlnw_dlnKA = dlnw_dlnKA(:);
dlnw_dlnKI = dlnw_dlnKI(:);

% Assemble the matrices

dlnw_dx = sparse([r; r_A; r_I], [m; m_A; m_I], [dlnw_da; -dlnw_dlnKA; -dlnw_dlnKI], nr, nm);

dlnw_dq = zeros(nr, q_info.q.number);
dlnw_dq(:,index.Keq) = diag(0.5 * h .* share_sum) * q_info.M_qKeqind_to_qallKeq;
dlnw_dq(:,index.KV)  = eye(nr);
dlnw_dq(sub2ind([nr,q_info.q.number], r,   index.KM(:))) = - dlnw_da - 0.5 * h(r) .* n .* share_sum(r);
dlnw_dq(sub2ind([nr,q_info.q.number], r_A, index.KA(:))) = dlnw_dlnKA;
dlnw_dq(sub2ind([nr,q_info.q.number], r_I, index.KI(:))) = dlnw_dlnKI;
