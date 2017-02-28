# @Author: Jenkins Alec <alec>
# @Date:   2017-01-24T17:52:32-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-28T14:15:55-08:00


import numpy as np
import json

filespec = 'Msfixed'
savepath = '/Users/alec/UCSB/cofeb_analysis_data/ta/cal_parameters_'+filespec+'.json'

# in CGS units
Mst = (6.54577e-5)
MstError = 1.01402555e-6
t = 1.0e-7
Ms = Mst/t
Keff = 8.20955e4
DW_width = 36.991e-7
Aex = 1.1233e-6

# radians and cm
height = 63.4248e-7 # calibrated height using +/- 3 points from cal peak and fitted finite edge width
heightError = 4.861e-7
theta = 1.033522705227361 # extracted from Bz0 and NV splitting away from edge
thetaError =  0.1089657# std of Bz0 and NV splitting away from edge
phi = 0*np.pi/180

# cm
scanSize = 5*(0.6e-4)

calparams = {'Ms':Ms, 't':t, 'MstError':MstError, 'theta':theta,
             'thetaError':thetaError, 'phi':phi, 'height':height, 'heightError':heightError, 'Keff':Keff,
             'Aex':Aex, 'DW_width':DW_width, 'scanSize':scanSize}

with open(savepath, 'w') as fwrite:
    json.dump(calparams, fwrite)
