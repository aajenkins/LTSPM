# @Author: Jenkins Alec <alec>
# @Date:   2017-01-24T17:52:32-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-04-09T20:32:43-07:00

import numpy as np
import json

filespec = 'Msfixed'
savepath = '/Users/alec/UCSB/cofeb_analysis_data/irmn/cal_parameters_'+filespec+'.json'

# SI units
Mst = 1.07365e-3
MstError = 6.298e-5
t = 0.911e-9
Ms = Mst/t
Keff = 5.9994e4
DW_energy_density = 0.00360166258672
DMI = 1.02e-4
Aex = ((DW_energy_density + np.pi*DMI)**2)/(16*Keff)
DW_width = DW_energy_density/(4*Keff)
# DW_width = np.sqrt(Aex/Keff)

# radians and m
height = 111.605264497e-9
heightError = 4.50487423163e-9 # std of calibration curves
theta = 0.985437137418 # extracted from Bz0 and NV splitting away from edge
thetaError =  0.0342669099991# std of Bz0 and NV splitting away from edge
phi = 180*np.pi/180 #symmetry direction of Bnv image

# cm
scanSize = 0.5*5e-6

calparams = {'Ms':Ms, 't':t, 'MstError':MstError, 'theta':theta,
             'thetaError':thetaError, 'phi':phi, 'height':height, 'heightError':heightError, 'Keff':Keff,
             'Aex':Aex, 'DW_width':DW_width, 'scanSize':scanSize}

with open(savepath, 'w') as fwrite:
    json.dump(calparams, fwrite)
