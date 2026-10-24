function d = cmb_basedir()

% d = cmb_basedir()
%
% Return the main directory of the model balancing repository

d = [fileparts(which(mfilename)) filesep '..' filesep '..' ];
