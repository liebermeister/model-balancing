% Model balancing: conversion between variables
%
%   mb_define_parameterisation - Define independent and dependent kinetic parameters
%   mb_kinetics_to_q           - Convert the kinetic constants in network.kinetics into a vector q of independent log kinetic constants
%   mb_qX_to_variables         - Compute all model variables from log kinetic constants q, log metabolite concentrations X, and fluxes V
%   mb_qX_to_y                 - Combine log kinetic constants q (vector) and log metabolite concentrations X (nm x ns matrix)
%   mb_q_to_kinetics           - Convert a vector q of independent log kinetic constants into a kinetics struct
%   mb_q_to_qall               - Convert a vector q of independent log kinetic constants into a vector qall of all log kinetic constants
%   mb_qall_to_kinetics        - Convert a vector qall of all log kinetic constants into a kinetics struct
%   mb_qall_to_q               - Not used (code is commented out)
%   mb_y_to_qX                 - Split a vector y of optimisation variables into log kinetic constants q (vector)
