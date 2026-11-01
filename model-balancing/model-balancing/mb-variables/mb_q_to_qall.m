function qall = mb_q_to_qall(q, q_info)

% qall = mb_q_to_qall(q, q_info)
%
% Convert a vector q of independent log kinetic constants into a vector qall of all log kinetic constants
% (entries of qall that depend on non-finite entries of q are set to nan)

ind_ok = isfinite(q);
qall   = q_info.M_q_to_qall(:,ind_ok) * q(ind_ok);
qall(sum(abs(q_info.M_q_to_qall(:,ind_ok)),2)==0) = nan;