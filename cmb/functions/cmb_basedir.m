function d = cmb_basedir()

d = [fileparts(which(mfilename)) filesep '..'  ];
