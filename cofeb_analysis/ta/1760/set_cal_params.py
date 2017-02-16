# @Author: Jenkins Alec <alec>
# @Date:   2017-01-24T17:52:32-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-14T21:47:27-08:00


import numpy as np

filespec = 'Msfixed'
savepath = '/Users/alec/UCSB/cofeb_analysis_data/ta/cal_parameters_'+filespec+'.txt'

# in CGS units
Mst = (6.22e-5)
t = 1.0e-7
Ms = Mst/t
theta = 55.58*np.pi/180
phi = 0*np.pi/180
h = 61.17 # in nm
herr = 4.33 # in nm
Keff = 8.82e4

calparams = [Ms, t, theta, phi, h, herr, Keff]

np.savetxt(savepath, calparams, delimiter=',')
