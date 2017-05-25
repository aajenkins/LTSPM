# @Author: Jenkins Alec <alec>
# @Date:   2017-01-24T17:52:32-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-04-10T09:48:40-07:00


import numpy as np
import json

filespec = 'Msfixed'
savepath = '/Users/alec/UCSB/cofeb_analysis_data/ta/cal_parameters_'+filespec+'.json'

# in Si units
Mst = (6.54577e-4)
MstError = 1.01402555e-5
t = 1.0e-9
Ms = Mst/t
Keff = 8.20955e3
DW_width = 36.991e-9
Aex = 1.1233e-11

# radians and m
height = 63.4248e-9 # calibrated height using +/- 3 points from cal peak and fitted finite edge width
heightError = 4.861e-9
theta = 1.033522705227361 # extracted from Bz0 and NV splitting away from edge
thetaError =  0.1089657# std of Bz0 and NV splitting away from edge
phi = 0

# m
scanSize = 0.6*(5e-6)

calparams = {'Ms':Ms, 't':t, 'MstError':MstError, 'theta':theta,
             'thetaError':thetaError, 'phi':phi, 'height':height, 'heightError':heightError, 'Keff':Keff,
             'Aex':Aex, 'DW_width':DW_width, 'scanSize':scanSize}

with open(savepath, 'w') as fwrite:
    json.dump(calparams, fwrite)
