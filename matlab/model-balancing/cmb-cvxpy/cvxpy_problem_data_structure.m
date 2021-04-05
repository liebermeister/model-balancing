function problem = cvxpy_problem_data_structure(network, q_info, prior, data, true_model)

% problem = cvxpy_problem_data_structure(network, q_info, prior, data, true_model)
%
% Convert a model balancing problem (network,q_info,true_model,prior,data) into a single data structure
% to be used by Elad's model balancing tool in python (based on CVXpy)
  
% Note:
% o activations and inhibitions are described by separate (positive) matrices of the shape #metabolites x #reactions
% o convention for order of KM, KA, and KI entries (in vectors) + their relation to entries in the matrix:
%   vector elements are ordered by their appearance in the respective (stoichiometric or regulation) matrix, in the format #reactions x # metabolites,
%   and following matlab's convention for indices in matrix elements (first column from top to bottom, then second column etc)
%   in python this corrresponds to matrix formatted as #metabolites x #reactions, and using the opposite convetion for indices 
%   (first row, then second row, etc)
% fields called "combined" contain the preposterior distribution of the respective quantities
  
% for the time being, data sets with changing flux directions cannot be handled:   

eval(default('true_model','[]'));  

if isempty(true_model),
  true_model.X = [];
  true_model.E = [];
  true_model.V = [];
  true_model.qall = nan * ones(q_info.qall.number,1);
end

M_q_to_qall = q_info.M_q_to_qall;

% If reactions need to be reoriented ..

if size(data.V.mean,2)>1,
if sum(max(sign(data.V.mean),2)-min(sign(data.V.mean),2)==2),
  error('data sets with changing flux directions cannot be handled');
else
  % reorient reactions and change all input data structures accordingly
  ind_neg = find(data.V.mean(:,1)<0);
  ind_pos = find(data.V.mean(:,1)>=0);
  
  % matrix for changes in qall parameter vector, due to reorientation:
  dum(q_info.qall.index.Keq) = 1;
  dum(q_info.qall.index.Keq(ind_neg)) = -1;
  dum(q_info.qall.index.KV) = 1;
  dum(q_info.qall.index.KM) = 1;
  dum(q_info.qall.index.KA) = 1;
  dum(q_info.qall.index.KI) = 1;
  dum(q_info.qall.index.Kcatf(ind_pos)) = 1;
  dum(q_info.qall.index.Kcatr(ind_pos)) = 1;
  Mreorient = diag(dum);
  Mreorient(q_info.qall.index.Kcatf(ind_neg),q_info.qall.index.Kcatr(ind_neg)) = eye(length(ind_neg));
  Mreorient(q_info.qall.index.Kcatr(ind_neg),q_info.qall.index.Kcatf(ind_neg)) = eye(length(ind_neg));

  M_q_to_qall = Mreorient * q_info.M_q_to_qall;
  
  % changes in network
  network.N(:,ind_neg)            = - network.N(:,ind_neg);
  network.kinetics.Keq(ind_neg)   = 1 ./ network.kinetics.Keq(ind_neg);
  dum                             = network.kinetics.Kcatf(ind_neg);
  network.kinetics.Kcatf(ind_neg) = network.kinetics.Kcatr(ind_neg); 
  network.kinetics.Kcatr(ind_neg) = dum; 

  % changes in true_model
  if ~isempty(true_model.V),
    true_model.V(ind_neg,:)               = - true_model.V(ind_neg,:);
    true_model.A_forward(ind_neg,:)       = - true_model.A_forward(ind_neg,:);
    true_model.kinetics.Keq(ind_neg)      = 1 ./ true_model.kinetics.Keq(ind_neg);
    dum                                   = true_model.kinetics.Kcatf(ind_neg);
    true_model.kinetics.Kcatf(ind_neg)    = true_model.kinetics.Kcatr(ind_neg); 
    true_model.kinetics.Kcatr(ind_neg)    = dum; 
  end
  
  % changes in data
  data.V.mean(ind_neg,:)    = - data.V.mean(ind_neg,:);
  data.qall.mean = Mreorient * data.qall.mean;
  
end
end

problem.standard_concentration        = '1 mM';
problem.network.metabolite_names      = network.metabolites;
problem.network.reaction_names        = network.actions;
problem.network.stoichiometric_matrix = full(network.N);
problem.network.activation_matrix     = full([network.regulation_matrix .* [network.regulation_matrix >0]]');
problem.network.inhibition_matrix     = full(-[network.regulation_matrix .* [network.regulation_matrix <0]]');

index                             = q_info.qall.index;
all_kinetic_constants.names       = q_info.qall.names;
all_kinetic_constants.true        = true_model.qall;
all_kinetic_constants.prior.mean  = M_q_to_qall * prior.q.mean;
all_kinetic_constants.prior.cov   = M_q_to_qall * inv(prior.q.cov_inv) * M_q_to_qall';
all_kinetic_constants.prior.std   = sqrt(diag(all_kinetic_constants.prior.cov));
all_kinetic_constants.data        = data.qall;
% PROBLEM: currently there is a prior only for INDEPENDENT Keq, not dependent ones; fix this!

Keq.unit          = 'depends on reaction stoichiometry';
Keq.true          = exp(all_kinetic_constants.true(index.Keq));
Keq.prior_ln.mean = all_kinetic_constants.prior.mean(index.Keq);
Keq.prior_ln.cov  = all_kinetic_constants.prior.cov(index.Keq,index.Keq);
Keq.data_ln.mean  = all_kinetic_constants.data.mean(index.Keq);
Keq.data_ln.std   = all_kinetic_constants.data.std(index.Keq);
[Keq.combined.mean_ln, Keq.combined.cov_ln] = cvxpy_combine_prior_and_likelihood(Keq.prior_ln.mean,Keq.prior_ln.cov,Keq.data_ln.mean,Keq.data_ln.std);
Keq.combined.geom_mean = exp(Keq.combined.mean_ln);

Kcatf.unit          = '1/s';
Kcatf.true          = exp(all_kinetic_constants.true(index.Kcatf));
Kcatf.prior_ln.mean = all_kinetic_constants.prior.mean(index.Kcatf);
Kcatf.prior_ln.cov  = all_kinetic_constants.prior.cov(index.Kcatf,index.Kcatf);
Kcatf.data_ln.mean  = all_kinetic_constants.data.mean(index.Kcatf);
Kcatf.data_ln.std   = all_kinetic_constants.data.std(index.Kcatf);
[Kcatf.combined.mean_ln, Kcatf.combined.cov_ln] = cvxpy_combine_prior_and_likelihood(Kcatf.prior_ln.mean,Kcatf.prior_ln.cov,Kcatf.data_ln.mean,Kcatf.data_ln.std);
Kcatf.combined.geom_mean = exp(Kcatf.combined.mean_ln);

Kcatr.unit          = '1/s';
Kcatr.true          = exp(all_kinetic_constants.true(index.Kcatr));
Kcatr.prior_ln.mean = all_kinetic_constants.prior.mean(index.Kcatr);
Kcatr.prior_ln.cov  = all_kinetic_constants.prior.cov(index.Kcatr,index.Kcatr);
Kcatr.data_ln.mean  = all_kinetic_constants.data.mean(index.Kcatr);
Kcatr.data_ln.std   = all_kinetic_constants.data.std(index.Kcatr);
[Kcatr.combined.mean_ln, Kcatr.combined.cov_ln] = cvxpy_combine_prior_and_likelihood(Kcatr.prior_ln.mean,Kcatr.prior_ln.cov,Kcatr.data_ln.mean,Kcatr.data_ln.std);
Kcatr.combined.geom_mean = exp(Kcatr.combined.mean_ln);

KM.unit          = 'mM';
KM.true          = exp(all_kinetic_constants.true(index.KM));
KM.prior_ln.mean = all_kinetic_constants.prior.mean(index.KM);
KM.prior_ln.cov  = all_kinetic_constants.prior.cov(index.KM,index.KM);
KM.data_ln.mean  = all_kinetic_constants.data.mean(index.KM);
KM.data_ln.std   = all_kinetic_constants.data.std(index.KM);
[KM.combined.mean_ln, KM.combined.cov_ln] = cvxpy_combine_prior_and_likelihood(KM.prior_ln.mean,KM.prior_ln.cov,KM.data_ln.mean,KM.data_ln.std);
KM.combined.geom_mean = exp(KM.combined.mean_ln);

KA.unit          = 'mM';
KA.true          = exp(all_kinetic_constants.true(index.KA));
KA.prior_ln.mean = all_kinetic_constants.prior.mean(index.KA);
KA.prior_ln.cov  = all_kinetic_constants.prior.cov(index.KA,index.KA);
KA.data_ln.mean  = all_kinetic_constants.data.mean(index.KA);
KA.data_ln.std   = all_kinetic_constants.data.std(index.KA);
[KA.combined.mean_ln, KA.combined.cov_ln] = cvxpy_combine_prior_and_likelihood(KA.prior_ln.mean,KA.prior_ln.cov,KA.data_ln.mean,KA.data_ln.std);
KA.combined.geom_mean = exp(KA.combined.mean_ln);

KI.unit          = 'mM';
KI.true          = exp(all_kinetic_constants.true(index.KI));
KI.prior_ln.mean = all_kinetic_constants.prior.mean(index.KI);
KI.prior_ln.cov  = all_kinetic_constants.prior.cov(index.KI,index.KI);
KI.data_ln.mean  = all_kinetic_constants.data.mean(index.KI);
KI.data_ln.std   = all_kinetic_constants.data.std(index.KI);
[KI.combined.mean_ln, KI.combined.cov_ln] = cvxpy_combine_prior_and_likelihood(KI.prior_ln.mean,KI.prior_ln.cov,KI.data_ln.mean,KI.data_ln.std);
KI.combined.geom_mean = exp(KI.combined.mean_ln);

problem.kinetic_constants.Keq   = Keq;
problem.kinetic_constants.Kcatf = Kcatf;
problem.kinetic_constants.Kcatr = Kcatr;
problem.kinetic_constants.KM    = KM;
problem.kinetic_constants.KA    = KA;
problem.kinetic_constants.KI    = KI;


% --- Metabolite concentrations

[nm,nc] = size(data.X.mean);

met.unit     = 'mM'; 
met.true     = exp(true_model.X);
met.prior_ln = prior.X; 
met.data_ln  = data.X;

% combine prior and data
[met_mean_ln, met_cov_ln] = cvxpy_combine_prior_and_likelihood(met.prior_ln.mean(:),diag(met.prior_ln.std(:).^2),met.data_ln.mean(:),met.data_ln.std(:));
met.combined.geom_mean  = exp(reshape(met_mean_ln,nm,nc));
met.combined.geom_std   = exp(reshape(sqrt(diag(met_cov_ln)),nm,nc));

problem.metabolite_concentrations = met;


% --- Enzyme concentrations

[nr,nc] = size(data.lnE.mean);

enz.unit         = 'mM';
enz.true         = true_model.E;
enz.prior_ln     = prior.lnE;

% convert (non-logarithmic) enzyme arithmetic mean + std (for prior and data) to logarithmic mean + std
% [enz_prior_mean_ln,enz_prior_std_ln] = lognormal_normal2log(enz.prior.mean(:),enz.prior.std(:));
% enz.data_ln      = data.lnE;
% [enz_data_mean_ln,enz_data_std_ln] = lognormal_normal2log(enz.data.mean(:),enz.data.std(:));

enz_prior_mean_ln = enz.prior_ln.mean(:);
enz_prior_std_ln  = enz.prior_ln.std(:);
enz.data_ln       = data.lnE;
enz_data_mean_ln  = enz.data_ln.mean(:);
enz_data_std_ln   = enz.data_ln.std(:);

% combine prior and data
[enz_mean_ln, enz_cov_ln] = cvxpy_combine_prior_and_likelihood(enz_prior_mean_ln,diag(enz_prior_std_ln.^2),enz_data_mean_ln,enz_data_std_ln);
enz.combined.geom_mean = exp(reshape(enz_mean_ln,nr,nc));
enz.combined.geom_std  = exp(reshape(sqrt(diag(enz_cov_ln)),nr,nc));

problem.enzyme_concentrations = enz;

% --- Reaction fluxes

problem.reaction_fluxes.unit      = 'mM/s';
problem.reaction_fluxes.true      = true_model.V;
problem.reaction_fluxes.data      = data.V;
