function pp = cmb_make_pp(network)

% Structure pp contains the network and other information (needed to call ecm enzyme cost score function)

pp.network             = network;
pp.enzyme_cost_weights = ones(size(network.actions));
pp.ind_scored_enzymes  = [1:length(network.actions)]';
pp.multiple_conditions = 0; 
pp.fluctuations_safety_margin = 0;
