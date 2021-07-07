function problem = cvxpy_problem_data_structure(network, q_info, prior, data, bounds, true_model, cmb_options)

% problem = cvxpy_problem_data_structure(network, q_info, prior, data, bounds, true_model, cmb_options)
%
% Convert a model balancing problem (data structures: network, q_info, true_model, prior, data) 
% into a single data structure, exported to json and then used by python model balancing tool (based on CVXpy)
%
% Structure, default values, and matrix sizes for a network with nr reactions and nm metabolites (and ns metabolic states)
%   .standard_concentration: '1 mM'
%   .state_names:            {ns x 1 cell}
%   .use_pseudo_values:      Boolean
%
%   .network
%   .network.metabolite_names:      {nm x 1 cell}
%   .network.reaction_names:        {nr x 1 cell}
%   .network.stoichiometric_matrix: [nm x nr double]  (sparse matrix)
%   .network.activation_matrix:     [nm x nr double]  (sparse matrix)
%   .network.inhibition_matrix:     [nm x nr double]  (sparse matrix)
%
%   .kinetic_constants
%   .kinetic_constants.Keq
%   .kinetic_constants.Kcatf
%   .kinetic_constants.Kcatr
%   .kinetic_constants.KM
%   .kinetic_constants.KA
%   .kinetic_constants.KI
%
%     .. where each subfield of .kinetic_constants is a struct of the form
%       .unit:              string, typically 'mM' and '1/s'
%       .true:              (column vector or matrix; true values in the case of artificial data; otherwise empty)
%       .pprior_ln.mean     pseudoprior, mean vector of log values
%       .pprior_ln.cov      pseudoprior, covariance matrix of log values
%       .bounds_ln.min      lower bounds on log values
%       .bounds_ln.max      upper bounds on log values
%       .data_ln.mean       data mean vector of log values
%       .data_ln.cov        data covariance matrix of log values
%       .bounds.min         lower bounds on log values
%       .bounds.max         upper bounds on log values
%       .combined.geom_mean preposterior geometric mean vector
%       .combined.mean_ln   preposterior mean vector of log values
%       %.combined.cov_ln   preposterior covariance matrix of log values
%
%   .kinetic_constants.all
%   .kinetic_constants.all.names      parameter names, list of strings
%   .kinetic_constants.all.geom_mean  preposterior geometric mean vector
%   .kinetic_constants.all.mean_ln    preposterior mean vector of log values
%   .kinetic_constants.all.prec_ln    preposterior precision matrix of log values
%
%   .metabolite_concentrations
%   .metabolite_concentrations.unit:              string, default 'mM'
%   .metabolite_concentrations.true:              matrix of true values (only in models with artificial data; otherwise [])
%   .metabolite_concentrations.prior_ln.mean      matrix of prior mean values for metabolite log-concentrations
%   .metabolite_concentrations.prior_ln.std       matrix of prior std dev for metabolite log-concentrations
%   .metabolite_concentrations.bounds_ln.min      matrix of lower bounds on metabolite log-concentrations
%   .metabolite_concentrations.bounds_ln.max      matrix of upper bounds on metabolite log-concentrations
%   .metabolite_concentrations.data_ln.mean       matrix of data  mean values for metabolite log-concentrations
%   .metabolite_concentrations.data_ln.std        matrix of data  std dev for metabolite log-concentrations
%   .metabolite_concentrations.bounds.min         matrix of preposterior geom mean for metabolite concentrations
%   .metabolite_concentrations.bounds.max         matrix of preposterior geom std  for metabolite concentrations
%   .metabolite_concentrations.combined.geom_mean matrix of preposterior geom mean for metabolite concentrations
%   .metabolite_concentrations.combined.geom_std  matrix of preposterior geom std  for metabolite concentrations
%
%   .enzyme_concentrations
%    ...
%    (structure analogous to .metabolite_concentrations; note that lower bounds for log(e) are typically -inf)
%
%   .reaction_fluxes
%   .reaction_fluxes.unit:          string, default 'mM/s'
%   .reaction_fluxes.true:          (matrix of true values in the case of artificial data; otherwise empty)
%   .reaction_fluxes.data.mean      matrix of flux data mean values
%   .reaction_fluxes.data.std       matrix of flux data std values
%
% 
% Note:
%   o activations and inhibitions are described by separate (positive) matrices of the shape #metabolites x #reactions
%   o convention for order of KM, KA, and KI entries (in vectors) + their relation to entries in the matrix:
%     vector elements are ordered by their appearance in the respective (stoichiometric or regulation) matrix, in the format #reactions x # metabolites,
%     and following matlab's convention for indices in matrix elements (first column from top to bottom, then second column etc)
%     in python this corrresponds to matrix formatted as #metabolites x #reactions, and using the opposite convetion for indices 
%     (first row, then second row, etc)
%   fields called "combined" contain the preposterior distribution of the respective quantities
  
% for the time being, data sets with changing flux directions cannot be handled:   

eval(default('true_model','[]'));  

if isempty(true_model),
  true_model.X = [];
  true_model.E = [];
  true_model.V = [];
  true_model.qall = nan * ones(q_info.qall.number,1);
end

M_q_to_qall    = q_info.M_q_to_qall;
M_qdep_to_qall = q_info.M_qdep_to_qall;
index          = q_info.qall.index;

% If reactions need to be reoriented ..

if size(data.V.mean,2)>1,
if sum(max(sign(data.V.mean),[],2)-min(sign(data.V.mean),[],2)==2),
  warning('Data set with flux reversal: note that this cannot be handled by the convex pathyon solver');
else
  % reorient reactions and change all input data structures accordingly
  ind_neg = find(data.V.mean(:,1)<0);
  ind_pos = find(data.V.mean(:,1)>=0);
  
  % matrix for changes in qall parameter vector, due to reorientation:
  signs = [];
  signs(index.Keq) = 1;
  signs(index.Keq(ind_neg)) = -1;
  signs(index.KV) = 1;
  signs(index.KM) = 1;
  signs(index.KA) = 1;
  signs(index.KI) = 1;
  signs(index.Kcatf(ind_pos)) = 1;
  signs(index.Kcatr(ind_pos)) = 1;
  Mreorient = diag(signs);
  Mreorient(index.Kcatf(ind_neg),index.Kcatr(ind_neg)) = eye(length(ind_neg));
  Mreorient(index.Kcatr(ind_neg),index.Kcatf(ind_neg)) = eye(length(ind_neg));

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
problem.state_names                   = data.samples;
problem.use_pseudo_values             = cmb_options.use_pseudo_values;
problem.network.metabolite_names      = network.metabolites;
problem.network.reaction_names        = network.actions;
problem.network.stoichiometric_matrix = full(network.N);
problem.network.activation_matrix     = full([network.regulation_matrix .* [network.regulation_matrix >0]]');
problem.network.inhibition_matrix     = full(-[network.regulation_matrix .* [network.regulation_matrix <0]]');

%  ---------------------------------
% kinetic constants - combine prior (for basic variables) and pseudo values (for dependent variables)
% into a "pseudoprior" for all variables

qall_struct.names        = q_info.qall.names;
qall_struct.true         = true_model.qall;
% if cmb_options.use_pseudo_values,
%   qall_struct.pprior.prec  = M_q_to_qall * prior.q.prec * M_q_to_qall' + M_qdep_to_qall * prior.qdep_pseudo.prec * M_qdep_to_qall';
%   qall_struct.pprior.mean  = qall_struct.pprior.prec \ [M_q_to_qall * prior.q.prec * prior.q.mean + M_qdep_to_qall * prior.qdep_pseudo.prec * prior.qdep_pseudo.mean];
%   qall_struct.pprior.cov   = inv(qall_struct.pprior.prec);
% else
%   qall_struct.pprior.cov  = M_q_to_qall * inv(prior.q.prec) * M_q_to_qall';
%   qall_struct.pprior.mean = M_q_to_qall * prior.q.mean;
% end
% qall_struct.pprior.std   = sqrt(diag(qall_struct.pprior.cov));
qall_struct.bounds.min   = bounds.q_all_min;
qall_struct.bounds.max   = bounds.q_all_max;

% hide data that are not supposed to be known 
switch cmb_options.use_kinetic_data,
  case 'all',      remove_data = [];
  case 'only_Keq', remove_data = setdiff(1:length(data.qall.mean),index.Keq);
  case 'none',     remove_data = 1:length(data.qall.mean);
  otherwise, error('incorrect option');
end
qall_struct.data                   = data.qall;
qall_struct.data.mean(remove_data) = nan;
qall_struct.data.std(remove_data)  = nan;

% preposterior

%[qall_struct.combined.mean_ln, qall_struct.combined.cov_ln, ~, qall_struct.combined.prec_ln] = cvxpy_combine_prior_and_likelihood(qall_struct.pprior.mean,qall_struct.pprior.cov,qall_struct.data.mean,qall_struct.data.std);

preposterior = cmb_prepare_posterior(prior, data, cmb_options, q_info);

qall_struct.combined.mean_ln = preposterior.qall.mean;
qall_struct.combined.cov_ln  = preposterior.qall.cov;
qall_struct.combined.prec_ln = preposterior.qall.prec;

% --------------------------------------------------------------
% copy results into "problem" struct


problem.kinetic_constants.all.names     = qall_struct.names;
problem.kinetic_constants.all.mean_ln   = qall_struct.combined.mean_ln;
problem.kinetic_constants.all.prec_ln   = qall_struct.combined.prec_ln;
problem.kinetic_constants.all.geom_mean = exp(qall_struct.combined.mean_ln);

% PROBLEM: currently there is a prior only for INDEPENDENT Keq, not dependent ones; fix this!

Keq.unit               = 'depends on reaction stoichiometry';
Keq.true               = exp(qall_struct.true(index.Keq));
%Keq.pprior_ln.mean     = qall_struct.pprior.mean(index.Keq);
%Keq.pprior_ln.cov      = qall_struct.pprior.cov(index.Keq,index.Keq);
Keq.bounds_ln.min      = qall_struct.bounds.min(index.Keq);
Keq.bounds_ln.max      = qall_struct.bounds.max(index.Keq);
Keq.data_ln.mean       = qall_struct.data.mean(index.Keq);
Keq.data_ln.std        = qall_struct.data.std(index.Keq);
Keq.bounds.min         = exp(qall_struct.bounds.min(index.Keq));
Keq.bounds.max         = exp(qall_struct.bounds.max(index.Keq));
Keq.combined.mean_ln   = preposterior.qtypes.Keq.mean;
Keq.combined.prec_ln   = preposterior.qtypes.Keq.prec;
Keq.combined.geom_mean = exp(Keq.combined.mean_ln);

Kcatf.unit          = '1/s';
Kcatf.true          = exp(qall_struct.true(index.Kcatf));
%Kcatf.pprior_ln.mean = qall_struct.pprior.mean(index.Kcatf);
%Kcatf.pprior_ln.cov  = qall_struct.pprior.cov(index.Kcatf,index.Kcatf);
Kcatf.bounds_ln.min = qall_struct.bounds.min(index.Kcatf);
Kcatf.bounds_ln.max = qall_struct.bounds.max(index.Kcatf);
Kcatf.data_ln.mean  = qall_struct.data.mean(index.Kcatf);
Kcatf.data_ln.std   = qall_struct.data.std(index.Kcatf);
Kcatf.bounds.min    = exp(qall_struct.bounds.min(index.Kcatf));
Kcatf.bounds.max    = exp(qall_struct.bounds.max(index.Kcatf));
Kcatf.combined.mean_ln   = preposterior.qtypes.Kcatf.mean;
Kcatf.combined.prec_ln   = preposterior.qtypes.Kcatf.prec;
Kcatf.combined.geom_mean = exp(Kcatf.combined.mean_ln);

Kcatr.unit          = '1/s';
Kcatr.true          = exp(qall_struct.true(index.Kcatr));
%Kcatr.pprior_ln.mean = qall_struct.pprior.mean(index.Kcatr);
%Kcatr.pprior_ln.cov  = qall_struct.pprior.cov(index.Kcatr,index.Kcatr);
Kcatr.bounds_ln.min = qall_struct.bounds.min(index.Kcatr);
Kcatr.bounds_ln.max = qall_struct.bounds.max(index.Kcatr);
Kcatr.data_ln.mean  = qall_struct.data.mean(index.Kcatr);
Kcatr.data_ln.std   = qall_struct.data.std(index.Kcatr);
Kcatr.bounds.min    = exp(qall_struct.bounds.min(index.Kcatr));
Kcatr.bounds.max    = exp(qall_struct.bounds.max(index.Kcatr));
Kcatr.combined.mean_ln   = preposterior.qtypes.Kcatr.mean;
Kcatr.combined.prec_ln   = preposterior.qtypes.Kcatr.prec;
Kcatr.combined.geom_mean = exp(Kcatr.combined.mean_ln);

KM.unit          = 'mM';
KM.true          = exp(qall_struct.true(index.KM));
%KM.pprior_ln.mean = qall_struct.pprior.mean(index.KM);
%KM.pprior_ln.cov  = qall_struct.pprior.cov(index.KM,index.KM);
KM.bounds_ln.min = qall_struct.bounds.min(index.KM);
KM.bounds_ln.max = qall_struct.bounds.max(index.KM);
KM.data_ln.mean  = qall_struct.data.mean(index.KM);
KM.data_ln.std   = qall_struct.data.std(index.KM);
KM.bounds.min    = exp(qall_struct.bounds.min(index.KM));
KM.bounds.max    = exp(qall_struct.bounds.max(index.KM));
KM.combined.mean_ln   = preposterior.qtypes.KM.mean;
KM.combined.prec_ln   = preposterior.qtypes.KM.prec;
KM.combined.geom_mean = exp(KM.combined.mean_ln);

KA.unit          = 'mM';
KA.true          = exp(qall_struct.true(index.KA));
%KA.pprior_ln.mean = qall_struct.pprior.mean(index.KA);
%KA.pprior_ln.cov  = qall_struct.pprior.cov(index.KA,index.KA);
KA.bounds_ln.min = qall_struct.bounds.min(index.KA);
KA.bounds_ln.max = qall_struct.bounds.max(index.KA);
KA.data_ln.mean  = qall_struct.data.mean(index.KA);
KA.data_ln.std   = qall_struct.data.std(index.KA);
KA.bounds.min    = exp(qall_struct.bounds.min(index.KA));
KA.bounds.max    = exp(qall_struct.bounds.max(index.KA));
KA.combined.mean_ln   = preposterior.qtypes.KA.mean;
KA.combined.prec_ln   = preposterior.qtypes.KA.prec;
KA.combined.geom_mean = exp(KA.combined.mean_ln);

KI.unit          = 'mM';
KI.true          = exp(qall_struct.true(index.KI));
%KI.pprior_ln.mean = qall_struct.pprior.mean(index.KI);
%KI.pprior_ln.cov  = qall_struct.pprior.cov(index.KI,index.KI);
KI.bounds_ln.min = qall_struct.bounds.min(index.KI);
KI.bounds_ln.max = qall_struct.bounds.max(index.KI);
KI.data_ln.mean  = qall_struct.data.mean(index.KI);
KI.data_ln.std   = qall_struct.data.std(index.KI);
KI.bounds.min    = exp(qall_struct.bounds.min(index.KI));
KI.bounds.max    = exp(qall_struct.bounds.max(index.KI));
KI.combined.mean_ln   = preposterior.qtypes.KI.mean;
KI.combined.prec_ln   = preposterior.qtypes.KI.prec;
KI.combined.geom_mean = exp(KI.combined.mean_ln);

problem.kinetic_constants.Keq   = Keq;
problem.kinetic_constants.Kcatf = Kcatf;
problem.kinetic_constants.Kcatr = Kcatr;
problem.kinetic_constants.KM    = KM;
problem.kinetic_constants.KA    = KA;
problem.kinetic_constants.KI    = KI;

% problem.kinetic_constants.prec_blocks.prec_Keq_ln       = qall_struct.combined.prec_ln(index.Keq,index.Keq);
% problem.kinetic_constants.prec_blocks.prec_KM_ln        = qall_struct.combined.prec_ln(index.KM,index.KM);
% problem.kinetic_constants.prec_blocks.prec_KA_ln        = qall_struct.combined.prec_ln(index.KA,index.KA);
% problem.kinetic_constants.prec_blocks.prec_KI_ln        = qall_struct.combined.prec_ln(index.KI,index.KI);
% problem.kinetic_constants.prec_blocks.prec_Kcatf_ln     = qall_struct.combined.prec_ln(index.Kcatf,index.Kcatf);
% problem.kinetic_constants.prec_blocks.prec_Kcatr_ln     = qall_struct.combined.prec_ln(index.Kcatr,index.Kcatr);
% 
% problem.kinetic_constants.prec_blocks.prec_Keq_KM_ln    = qall_struct.combined.prec_ln(index.Keq,index.KM);
% problem.kinetic_constants.prec_blocks.prec_Keq_KA_ln    = qall_struct.combined.prec_ln(index.Keq,index.KA);
% problem.kinetic_constants.prec_blocks.prec_Keq_KI_ln    = qall_struct.combined.prec_ln(index.Keq,index.KI);
% problem.kinetic_constants.prec_blocks.prec_Keq_Kcatf_ln = qall_struct.combined.prec_ln(index.Keq,index.Kcatf);
% problem.kinetic_constants.prec_blocks.prec_Keq_Kcatr_ln = qall_struct.combined.prec_ln(index.Keq,index.Kcatr);
% 
% problem.kinetic_constants.prec_blocks.prec_KM_KA_ln    = qall_struct.combined.prec_ln(index.KM,index.KA);
% problem.kinetic_constants.prec_blocks.prec_KM_KI_ln    = qall_struct.combined.prec_ln(index.KM,index.KI);
% problem.kinetic_constants.prec_blocks.prec_KM_Kcatf_ln = qall_struct.combined.prec_ln(index.KM,index.Kcatf);
% problem.kinetic_constants.prec_blocks.prec_KM_Kcatr_ln = qall_struct.combined.prec_ln(index.KM,index.Kcatr);
% 
% problem.kinetic_constants.prec_blocks.prec_KA_KI_ln    = qall_struct.combined.prec_ln(index.KA,index.KI);
% problem.kinetic_constants.prec_blocks.prec_KA_Kcatf_ln = qall_struct.combined.prec_ln(index.KA,index.Kcatf);
% problem.kinetic_constants.prec_blocks.prec_KA_Kcatr_ln = qall_struct.combined.prec_ln(index.KA,index.Kcatr);
% 
% problem.kinetic_constants.prec_blocks.prec_KI_Kcatf_ln = qall_struct.combined.prec_ln(index.KI,index.Kcatf);
% problem.kinetic_constants.prec_blocks.prec_KI_Kcatr_ln = qall_struct.combined.prec_ln(index.KI,index.Kcatr);
% 
% problem.kinetic_constants.prec_blocks.prec_Kcatf_Kcatr_ln = qall_struct.combined.prec_ln(index.Kcatf,index.Kcatr);

% --- Metabolite concentrations

[nm,nc] = size(data.X.mean);

met.unit     = 'mM'; 
met.true     = exp(true_model.X);

met.bounds_ln.min = repmat(bounds.x_min,1,nc);
met.bounds_ln.max = repmat(bounds.x_max,1,nc);
met.prior_ln = prior.X;
met.data_ln  = data.X;
met.bounds.min         = exp(repmat(bounds.x_min,1,nc));
met.bounds.max         = exp(repmat(bounds.x_max,1,nc));
% combine prior and data
[met_mean_ln, met_cov_ln] = cvxpy_combine_prior_and_likelihood(met.prior_ln.mean(:),diag(met.prior_ln.std(:).^2),met.data_ln.mean(:),met.data_ln.std(:));
met.combined.geom_mean = exp(reshape(met_mean_ln,nm,nc));
met.combined.geom_std  = exp(reshape(sqrt(diag(met_cov_ln)),nm,nc));

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
enz.bounds_ln.min = repmat(log(bounds.e_min),1,nc);
enz.bounds_ln.max = repmat(log(bounds.e_max),1,nc);
enz.data_ln       = data.lnE;
enz_data_mean_ln  = enz.data_ln.mean(:);
enz_data_std_ln   = enz.data_ln.std(:);
enz.bounds.min         = exp(repmat(log(bounds.e_min),1,nc));
enz.bounds.max         = exp(repmat(log(bounds.e_max),1,nc));
% combine prior and data
[enz_mean_ln, enz_cov_ln] = cvxpy_combine_prior_and_likelihood(enz_prior_mean_ln,diag(enz_prior_std_ln.^2),enz_data_mean_ln,enz_data_std_ln);
enz.combined.geom_mean = exp(reshape(enz_mean_ln,nr,nc));
enz.combined.geom_std  = exp(reshape(sqrt(diag(enz_cov_ln)),nr,nc));

problem.enzyme_concentrations = enz;

% --- Reaction fluxes

problem.reaction_fluxes.unit      = 'mM/s';
problem.reaction_fluxes.true      = true_model.V;
problem.reaction_fluxes.data      = data.V;
