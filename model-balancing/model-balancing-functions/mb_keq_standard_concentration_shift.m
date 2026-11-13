function [shift_all, shift_ind] = mb_keq_standard_concentration_shift(network, q_info, mb_options)

% [shift_all, shift_ind] = mb_keq_standard_concentration_shift(network, q_info, mb_options)
%
% Shift of ln Keq that converts prior values and bounds (mb_options.quantities.Keq) from their standard
% concentration (mb_options.quantities.Keq.standard_concentration, by default the one of the prior table)
% to 1 mM, the standard concentration used internally in model balancing
%
% shift_all  one value per reaction (to be added to prior means and bounds of ln Keq)
% shift_ind  one value per independent equilibrium constant (vector q; mapped with q_info.M_qall_to_q)
%
% Water is not counted in the sums of stoichiometric coefficients (see standard_concentration_convert).
% Without the field quantities.Keq.standard_concentration, 1 mM is assumed (shift zero).

nr      = length(network.actions);
nKeqind = length(q_info.q.index.Keq);

shift_all = zeros(nr,1);
shift_ind = zeros(nKeqind,1);

if ~isfield(mb_options.quantities.Keq, 'standard_concentration')
  return
end

if standard_concentration_in_mM(mb_options.quantities.Keq.standard_concentration) == 1
  return
end

shift_all = standard_concentration_convert(zeros(nr,1), 'ln Keq', network, mb_options.quantities.Keq.standard_concentration, '1 mM');
shift_ind = full(q_info.M_qall_to_q(q_info.q.index.Keq, q_info.qall.index.Keq) * shift_all);
