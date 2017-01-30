import numpy as np

filespec = 'Msfixed'
savepath = '/Users/alec/UCSB/LTSPM/cofeb_analysis/irmn/cal_parameters_'+filespec+'.txt'

Ms = 1.486e6
t = 0.911
theta = 54.30*np.pi/180
phi = 180*np.pi/180
h = 141.06
herr = 6.25

calparams = [Ms, t, theta, phi, h, herr]

np.savetxt(savepath, calparams, delimiter=',')

filespec = 'Msnotfixed'
savepath = '/Users/alec/UCSB/LTSPM/cofeb_analysis/irmn/cal_parameters_'+filespec+'.txt'

Ms = 1.068e6
t = 1
theta = 55.01*np.pi/180
phi = 180*np.pi/180
h = 110.97
herr = 4.15

calparams = [Ms, t, theta, phi, h, herr]

np.savetxt(savepath, calparams, delimiter=',')
