function [q,X] = cmb_y_to_qX(y,nm,ns)

% [q, X] = cmb_y_to_qX(y, nm, ns)
%
% Split a vector y of optimisation variables into log kinetic constants q (vector)
% and log metabolite concentrations X (nm x ns matrix) (see also cmb_qX_to_y)

X = reshape(y(1:nm*ns),nm,ns);
q = y(nm*ns+1:end);
