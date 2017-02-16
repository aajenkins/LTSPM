# @Author: Jenkins Alec <alec>
# @Date:   2017-01-24T17:52:32-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-14T15:30:01-08:00


import numpy as np

filespec = 'Msfixed'
savepath = '/Users/alec/UCSB/cofeb_analysis_data/irmn/cal_parameters_'+filespec+'.txt'

Mst = (1.074e-4)
t = 0.911e-7
Ms = Mst/t
theta = 54.99*np.pi/180
phi = 180*np.pi/180
h = 111.6
herr = 4.5
Keff = 5.999e5

calparams = [Ms, t, theta, phi, h, herr, Keff]

np.savetxt(savepath, calparams, delimiter=',')
