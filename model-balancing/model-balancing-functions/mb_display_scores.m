function improvement = mb_display_scores(network, q_info, mb_options, pp, preposterior, init, optimal, true_model, V, verbose)

% improvement = mb_display_scores(network, q_info, mb_options, pp, preposterior, init, optimal, true_model, V, verbose)
%
% Display the log posterior scores for the true values (if given), the initial values, and the optimised values
%
% Output variables
%   improvement  Boolean, true if the optimised values have a higher log posterior than the initial values

[nm,nr] = size(network.N);

eval(default('true_model','[]','verbose','0'));

if ~isempty(true_model),
  y_true     = mb_qX_to_y(true_model.q,true_model.X,nm,mb_options.ns);
  score_true = mb_log_posterior(y_true, pp, preposterior, V, mb_options, q_info, verbose);
  display(sprintf('log posterior TRUE:      %f', score_true ));
end

y_init     = mb_qX_to_y(init.q, init.X,nm,mb_options.ns);
score_init = mb_log_posterior(y_init, pp, preposterior, V, mb_options, q_info, verbose);
display(sprintf('log posterior INITIAL:   %f', score_init ));

y_opt     = mb_qX_to_y(optimal.q,optimal.X,nm,mb_options.ns);
score_opt =  mb_log_posterior(y_opt, pp, preposterior, V, mb_options, q_info, verbose);
display(sprintf('log posterior OPTIMISED: %f', score_opt ));

% did the last run improve the objective? (assuming that the scores are negative)
if score_opt<0,
  improvement = [score_opt > 0.999999999 * score_init];
else
  improvement = [score_opt > 1.000000001 * score_init];
end