function y = cmb_qX_to_y(q,X,nm,ns)

y(1:nm*ns,1) = reshape(X,nm*ns,1);
y(nm*ns+[1:length(q)],1) = q;
