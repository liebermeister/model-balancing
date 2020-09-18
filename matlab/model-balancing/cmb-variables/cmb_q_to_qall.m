function qall = cmb_q_to_qall(q, q_info)

ind_ok = isfinite(q);
qall   = q_info.M_q_to_qall(:,ind_ok) * q(ind_ok);
qall(sum(abs(q_info.M_q_to_qall(:,ind_ok)),2)==0) = nan;