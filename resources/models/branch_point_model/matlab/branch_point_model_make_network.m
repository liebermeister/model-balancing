% -------------------------------------------------------------
% Branch point model (network structure from SBML or SBtab file)
%
% X1       
%   \v1 
%    \    v3
%      S4 --> X3
%    /
%   /v2
% X2
%
% -------------------------------------------------------------

% -------------------------------------------------------------
% Define the network model

metabolites = {'X1','X2','X3','S4'}';

N = [-1  0  0;...
      0 -1  0 ;...
      0  0  1;...
      1  1 -1];

reversible = ones(3,1);
external   = [1 2 3]';
network    = network_construct(N,reversible,external,metabolites);

sbml_file = [cmb_resourcedir '/models/branch_point_model/branch_point_model.xml'];
network_sbml_export(network, 0, 'branch_point_model', sbml_file);

% graphics

layout_file = [cmb_resourcedir '/models/branch_point_model/branch_point_model_Position.tsv'];
network = netgraph_read_positions(network, layout_file);
network.graphics_par.squaresize=0.02;

% figure(1000)
% network = netgraph_edit_positions(network,layout_file); 

figure(1000)
netgraph_concentrations(network,network.external,[],1);
graphics_file = [cmb_resourcedir '/models/branch_point_model/graphics/branch_point_model.eps'];
print(graphics_file, '-f1000', '-depsc'); 

% export sbtab model file
sbtab_file = [cmb_resourcedir '/models/branch_point_model/branch_point_model.tsv'];
network_to_sbtab(network, struct('filename', sbtab_file));
