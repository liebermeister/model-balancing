function improvement = cmb_display_scores(network, q_info, cmb_options, pp, preposterior, init, optimal, true, V, verbose)

[nm,nr] = size(network.N);

eval(default('true','[]','verbose','0'));

if ~isempty(true),
  y_true     = cmb_qX_to_y(true.q,true.X,nm,cmb_options.ns);
  score_true = cmb_log_posterior(y_true,pp,preposterior,true.V,cmb_options,q_info, verbose);
  display(sprintf('log posterior TRUE:      %f', score_true ));
end

y_init     = cmb_qX_to_y(init.q, init.X,nm,cmb_options.ns);
score_init = cmb_log_posterior(y_init, pp, preposterior, V, cmb_options, q_info, verbose);
display(sprintf('log posterior INITIAL:   %f', score_init ));

y_opt     = cmb_qX_to_y(optimal.q,optimal.X,nm,cmb_options.ns);
score_opt =  cmb_log_posterior(y_opt, pp, preposterior, V, cmb_options, q_info, verbose);
display(sprintf('log posterior OPTIMISED: %f', score_opt ));

% did the last run improve the objective? (assuming that the scores are negative)
if score_opt<0,
  improvement = [score_opt > 0.999999999 * score_init];
else
  improvement = [score_opt > 1.000000001 * score_init];
end