load ./funcs/mat_files/Cdisk
load ./funcs/mat_files/params

% verify if LPVcore is in the path (via pmatrix object)
assert(exist('pmatrix','class') == 8, ...
    'LPVcore not (fully) included in the path')

addpath('funcs/Simulator/Base Simulator')
addpath('funcs/misc/')

% Load Gyroscope parameters
load ./funcs/Simulator/params

% Physical parameters
Ib = params.Ib;     Ic = params.Ic;     Id = params.Id;
Jb = params.Jb;     Jc = params.Jc;     Jd = params.Jd;
Ka = params.Ka;     Kb = params.Kb;     Kc = params.Kc;
Km1 = params.Km1;   Km2 = params.Km2;
Km3 = params.Km3;   Km4 = params.Km4;
fv1 = params.fv1;   fv2 = params.fv2;
fv3 = params.fv3;   fv4 = params.fv4;

% Reference low pass filter parameters
lowpass_num = [1];
lowpass_den = [1];
