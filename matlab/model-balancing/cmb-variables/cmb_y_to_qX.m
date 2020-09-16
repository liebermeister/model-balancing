function [q,X] = cmb_y_to_qX(y,nm,ns)

X = reshape(y(1:nm*ns),nm,ns);
q = y(nm*ns+1:end);
