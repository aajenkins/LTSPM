import numpy as np

filenum = 1760

filespec = 'Msfixed'
savepath = '/Users/alec/UCSB/LTSPM/cofeb_analysis/ta/'+str(filenum)+'/cal_parameters_'+filespec+'.txt'

Ms = 6.496e5
t = 1
theta = 55.61*np.pi/180
phi = 0*np.pi/180
h = 63.33
herr = 4.53

calparams = [Ms, t, theta, phi, h, herr]

np.savetxt(savepath, calparams, delimiter=',')

filespec = 'Msnotfixed'
savepath = '/Users/alec/UCSB/LTSPM/cofeb_analysis/ta/'+str(filenum)+'/cal_parameters_'+filespec+'.txt'

Ms = 7.72e5
t = 1
theta = 55.72*np.pi/180
phi = 0*np.pi/180
h = 72.34
herr = 2.79

calparams = [Ms, t, theta, phi, h, herr]

np.savetxt(savepath, calparams, delimiter=',')