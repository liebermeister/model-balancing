% -------------------------------------------------------------
% Model and data for double branch example model
%
% X1 <-> S2 <-> S3 <-> X4
% -------------------------------------------------------------

% -------------------------------------------------------------
% Define the network model
  
metabolites = {'X1','S2','X3','X4'}';

N = [ -1  0  0; ...
       1 -1  0; ...
       0  1 -1; ...
       0  0  1]; 

reversible = [1 1 1]';
external   = [1 4]';
network    = network_construct(N,reversible,external,metabolites);

sbml_file = [cmb_resourcedir '/models/three_chain_model/three_chain_model.xml'];
network_sbml_export(network, 0, 'three_chain_model', sbml_file);

% graphics

layout_file = [cmb_resourcedir '/models/three_chain_model/three_chain_model_Position.tsv'];
network = netgraph_read_positions(network, layout_file);
network.graphics_par.squaresize=0.01;

% figure(1000)
% network = netgraph_edit_positions(network,layout_file); 

figure(1000)
netgraph_concentrations(network,network.external,[],1);
graphics_file = [cmb_resourcedir '/models/three_chain_model/graphics/three_chain_model.eps'];
print(graphics_file, '-f1000', '-depsc'); 

% export sbtab model file 
sbtab_file = [cmb_resourcedir '/models/three_chain_model/three_chain_model.tsv'];
network_to_sbtab(network, struct('filename', sbtab_file));

