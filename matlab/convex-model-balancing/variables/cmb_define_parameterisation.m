function q_info = cmb_define_parameterisation(network, cmb_options)

% q_info = cmb_define_parameterisation(network, cmb_options)
%
% q_info is struct containing a description of the indepdendent / all parameters
%                  as well as the conversion matrices between them 
%
% Parameterisation 'KV_KM_KA_KI':
%  q contains independent parameters: KV, KM, KA, KI
%  qall contains all parameters:      KV, KM, KA, KI
%
% Parameterisation 'Keq_KV_KM_KA_KI':
%  q contains independent parameters: Keq_ind, KV, KM, KA, KI
%  qall contains all parameters:      Keq,     KV, KM, KA, KI, Kcat_f, Kcat_r
%
% Defining the independent equilibrium constants:  
%   split Ntot = A * B, where A has full column rank
%   then we know that we can write keq = B' * keqind, where keqind contains independent constants
  
[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);

h = ones(nr,1);

switch cmb_options.parameterisation,  
  
  case 'KV_KM_KA_KI', 
    %% assuming Keq are fixed, as given by model

    %% Independent parameters
    
     q.names = [repmat({'KV'},nr,1); ...
                repmat({'KM'},nKM,1); ...
                repmat({'KA'},nKA,1); ...
                repmat({'KI'},nKI,1) ];
     
     q.index.KV  = 1:nr;    
     q.index.KM  = nr + [1:nKM];
     q.index.KA  = nr + nKM + [1:nKA];
     q.index.KI  = nr + nKM + nKA + [1:nKI];

     %% All parameters
     qall = q;
     
  case 'Keq_KV_KM_KA_KI',
    %% include independent Keq values into parameter vector 
    
    %%%% --------- NOT YET SURE IF THIS IS CORRECT .. CHECK AGAIN !!!!
    [Bt, At, independent_reactions] = reduce_N(network.N');
    nKeqind = size(Bt,2);
    
    %% Independent parameters
    
    q.names = [repmat({'Keqind'},nKeqind,1); ...
               repmat({'KV'},nr,1); ...
               repmat({'KM'},nKM,1); ...
               repmat({'KA'},nKA,1); ...
               repmat({'KI'},nKI,1) ];
    
    q.index.Keq = 1:nKeqind;
    q.index.KV  = nKeqind + [1:nr];    
    q.index.KM  = nKeqind + nr + [1:nKM];
    q.index.KA  = nKeqind + nr + nKM + [1:nKA];
    q.index.KI  = nKeqind + nr + nKM + nKA + [1:nKI];
    
     %% All parameters

     qall.names = [repmat({'Keq'},nr,1); ...
                   repmat({'KV'},nr,1); ...
                   repmat({'KM'},nKM,1); ...
                   repmat({'KA'},nKA,1); ...
                   repmat({'KI'},nKI,1);
                   repmat({'Kcatf'},nr,1); ...
                   repmat({'Kcatr'},nr,1); ...
]; 
     
     qall.index.Keq   = 1:nr;
     qall.index.KV    = nr + [1:nr];    
     qall.index.KM    = 2 * nr + [1:nKM];
     qall.index.KA    = 2 * nr + nKM + [1:nKA];
     qall.index.KI    = 2 * nr + nKM + nKA + [1:nKI];
     qall.index.Kcatf = 2 * nr + nKM + nKA + nKI + [1:nr];    
     qall.index.Kcatr = 2 * nr + nKM + nKA + nKI + nr + [1:nr];    
     
  otherwise error('Parameterisation not supported');

end

q.number     = length(q.names);
qall.number = length(qall.names);
     


% ---------------------------------------------------
% Conversion matrices

M_q_to_qall = zeros(qall.number, q.number);
M_qall_to_q = zeros(q.number, qall.number);

M_q_to_qall(qall.index.KV, q.index.KV) = eye(nr);
M_q_to_qall(qall.index.KM, q.index.KM) = eye(nKM);
M_q_to_qall(qall.index.KA, q.index.KA) = eye(nKA);
M_q_to_qall(qall.index.KI, q.index.KI) = eye(nKI);

M_qall_to_q(q.index.KV, qall.index.KV) = eye(nr);
M_qall_to_q(q.index.KM, qall.index.KM) = eye(nKM);
M_qall_to_q(q.index.KA, qall.index.KA) = eye(nKA);
M_qall_to_q(q.index.KI, qall.index.KI) = eye(nKI);

switch cmb_options.parameterisation,  
  case 'Keq_KV_KM_KA_KI',
    %% formula for kcat+- (see Eq. 25 in modular rate law paper)
    stoich_matrix_reaction_KM = zeros(nr,nKM);
    [ind_r,ind_m] = ind2sub([nr,nm],KM_ind);
    for it = 1:nKM,
      stoich_matrix_reaction_KM(ind_r(it),it) = network.N(ind_m(it),ind_r(it));
    end
    M_q_to_qall(qall.index.Keq, q.index.Keq)  = Bt;
    M_q_to_qall(qall.index.Kcatf, q.index.Keq)=  0.5 * diag(h) * Bt;
    M_q_to_qall(qall.index.Kcatr, q.index.Keq)= -0.5 * diag(h) * Bt;
    M_q_to_qall(qall.index.Kcatf, q.index.KM) =  0.5 * diag(h) * stoich_matrix_reaction_KM;
    M_q_to_qall(qall.index.Kcatr, q.index.KM) = -0.5 * diag(h) * stoich_matrix_reaction_KM;
    M_q_to_qall(qall.index.Kcatf, q.index.KV) = eye(nr);
    M_q_to_qall(qall.index.Kcatr, q.index.KV) = eye(nr);

    M_qall_to_q(q.index.Keq, qall.index.Keq)  = pinv(full(Bt));
    
end

q_info.q           = q;            % independent parameters
q_info.qall        = qall;         % all parameters
q_info.M_q_to_qall = M_q_to_qall;  % conversion independent -> all
q_info.M_qall_to_q = M_qall_to_q;  % conversion all -> independent

switch cmb_options.parameterisation,  
  case 'Keq_KV_KM_KA_KI',
    q_info.M_qallKeq_to_qKeqind = q_info.M_qall_to_q(q_info.q.index.Keq,q_info.qall.index.Keq);
    q_info.M_qKeqind_to_qallKeq = q_info.M_q_to_qall(q_info.qall.index.Keq,q_info.q.index.Keq);
end

if 0,
  % plot conversion matrices
  figure(1); im(q_info.M_qall_to_q,[],q_info.q.names,q_info.qall.names)   
  figure(2); im(q_info.M_q_to_qall,[],q_info.qall.names,q_info.q.names)
end

q_info.KM_matrix_indices = KM_ind;