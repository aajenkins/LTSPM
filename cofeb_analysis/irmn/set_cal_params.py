# @Author: Jenkins Alec <alec>
# @Date:   2017-01-24T17:52:32-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-21T19:31:13-08:00


import numpy as np

filespec = 'Msfixed'
savepath = '/Users/alec/UCSB/cofeb_analysis_data/irmn/cal_parameters_'+filespec+'.txt'

# cgs units
Mst = 1.07365e-4
t = 0.911e-7
Ms = Mst/t
Keff = 5.9994e5
DW_width = 1.523e-6
Aex = 1.3919e-6

# radians and cm
theta = 54.992*np.pi/180
phi = 180*np.pi/180
h = 111.61e-7
herr = 4.5e-7

# cm
scanSize = 5*0.5e-4



calparams = [Ms, t, theta, phi, h, herr, Keff, Aex, DW_width, scanSize]

np.savetxt(savepath, calparams, delimiter=',')
