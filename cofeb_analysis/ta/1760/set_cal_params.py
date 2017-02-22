# @Author: Jenkins Alec <alec>
# @Date:   2017-01-24T17:52:32-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-21T20:08:37-08:00


import numpy as np

filespec = 'Msfixed'
savepath = '/Users/alec/UCSB/cofeb_analysis_data/ta/cal_parameters_'+filespec+'.txt'

# in CGS units
Mst = (6.54577e-5)
t = 1.0e-7
Ms = Mst/t
Keff = 8.20955e4
DW_width = 36.991e-7
Aex = 1.1233e-6

# radians and cm
h = 55.865e-7
herr = 4.632e-7
theta = 55.357*np.pi/180
phi = 0*np.pi/180

# cm
scanSize = 5*(0.6e-4)

calparams = [Ms, t, theta, phi, h, herr, Keff, Aex, DW_width, scanSize]

np.savetxt(savepath, calparams, delimiter=',')
