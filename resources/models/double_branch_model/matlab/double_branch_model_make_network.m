% -------------------------------------------------------------
% Model and data for double branch example model
%
%     E1
% X1  -> S2           -> X6
%          \E3    E4/
% 	      S5 
%     E2   /      E5\
% X3  -> S4           -> X7
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

sbml_file = [cmb_resourcedir '/models/double_branch_model/double_branch_model.xml'];
network_sbml_export(network, 0, 'double_branch_model', sbml_file);

% graphics
layout_file = [cmb_resourcedir '/models/double_branch_model/double_branch_model_Position.tsv'];
network = netgraph_read_positions(network, layout_file);
network.graphics_par.squaresize=0.05;

% figure(1000)
% network = netgraph_edit_positions(network,layout_file); 

figure(1000)
netgraph_concentrations(network,network.external,[],1);
graphics_file = [cmb_resourcedir '/models/double_branch_model/graphics/double_branch_model.eps'];
print(graphics_file, '-f1000', '-depsc'); 

% export sbtab model file 
sbtab_file = [cmb_resourcedir '/models/double_branch_model/double_branch_model.tsv'];
network_to_sbtab(network, struct('filename', sbtab_file));

