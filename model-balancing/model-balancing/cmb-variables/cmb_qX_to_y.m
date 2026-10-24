function y = cmb_qX_to_y(q,X,nm,ns)

% y = cmb_qX_to_y(q, X, nm, ns)
%
% Combine log kinetic constants q (vector) and log metabolite concentrations X (nm x ns matrix)
% into a single vector y of optimisation variables (see also cmb_y_to_qX)

y(1:nm*ns,1) = reshape(X,nm*ns,1);
y(nm*ns+[1:length(q)],1) = q;
