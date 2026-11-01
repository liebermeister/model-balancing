function q = mb_qall_to_q(qall,q_info)

% q = mb_qall_to_q(qall, q_info)
%
% Convert a vector qall of all log kinetic constants into a vector q of independent log kinetic constants
  
q = q_info.M_qall_to_q * qall;