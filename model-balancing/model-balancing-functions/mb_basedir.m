function d = mb_basedir()

% d = mb_basedir()
%
% Return the main directory of the model balancing repository

d = [fileparts(which(mfilename)) filesep '..' filesep '..' ];
