# ID_CONTROL_LPVSUBNET

Code associated with ...

## Required software

1.  LPVcore toolbox for MATLAB
1.  deepSI Python toolbox (https://github.com/GerbenBeintema/deepSI)
1.  The python scripts in 'autoencoders/' folders were executed using
Python 3.6.9 with keras 2.3.1 and tensorflow 1.14.0
1.  The python scripts in 'lpv_subnet/identification' folders were executed
using both
Python 3.8.10 with tensorflow 2.9.1 and Python 3.6.9 with tensorflow 2.6.21
1.  MATLAB version: R2022b/R2023a
1.  The MATLAB scripts in 'BLA/' folders were executed using the System
Identification Toolbox version: 9.14 (R2021a)

## Notes

1.  The code in 'lpv_subnet/' folders is based on
https://gitlab.com/releases-c-verhoek/lpv-subnet
1.  The code in 'autoencoders/' folders is based on
http://dysco.imtlucca.it/masti/autoencoders/

## Description

1.  'autoencoders/': identification using autoencoders; see
'gyroscope_main.py'.
1.  'BLA/': identification using BLA; see 'main.m'.
1.  'data/': folder containing data.
1.  'lpvcore/': identification using LPVcore; see 'gyroscope_main.m'.
1.  'lpvsubnet/identification': identification using LPVSUBNET;
see 'gyroscope.py';
use 'matrices.py' to export the LPVSUBNET to matlab
(creates a matlab file matrices.m that should be moved to 'lpvsubnet/control').
1.  'lpvsubnet/control': LPV controller design and simulation; see 'main.m'.
1.  'lpvsubnet/control_double_input': LPV controller design and simulation with
i_1 and i_2 as inputs; see 'main.m'. (NOT WORKING)
