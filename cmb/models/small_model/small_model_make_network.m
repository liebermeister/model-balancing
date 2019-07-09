% -------------------------------------------------------------
% Model and data for small example model
%
%     E1
% X1  -> S2           X6
%          \E3    E4/
% 	  --> S5 
%     E2   /      E5\
% X3  -> S4           X7
%
% -------------------------------------------------------------

% -------------------------------------------------------------
% Define the network model
  
metabolites = {'X1','S2','X3','S4','S5','X6','X7'}';

N = [-1 1 0 0 0 0 0 ;...
     0 0 -1 1 0 0 0 ;...
     0 -1 0 -1 1 0 0 ;...
     0 0 0 0 -1 1 0 ;
     0 0 0 0 -1 0 1]';

reversible = ones(5,1);
external   = [1 3 6 7]';
network    = network_construct(N,reversible,external,metabolites);

network_sbml_export(network, 0, 'small_model', [cmb_basedir '/models/small_model/small_model.xml']);

% graphics

figure(1000)
layout_file = [cmb_basedir '/models/small_model/small_model_Layout.tsv]'
network = netgraph_edit_positions(network,layout_file);%,flag_KEGG_ids,flag_element_names,flag_fixed,flag_save_metabolites_as_fixed);

figure(1000)
netgraph_concentrations(network,network.external,[],1);
print([cmb_basedir '/models/small_model/small_model.eps'], '-f1000', '-depsc'); 
