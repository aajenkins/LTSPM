# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:17:32 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import cwesr_fit_single as cwesr

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)

#material parameters
pi = np.pi

#file constants

basepath = '/Users/alec/UCSB/calibration_data/tacofeb_nv_angle/anglecal'
num_avg = 10
filenum = 2

data = np.loadtxt(basepath+str(filenum).zfill(6)+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]
    
edata = np.transpose(data)
#print(edata)
x, y = edata

cwresult = cwesr.cwesr_fit(x,y,gauss=False,gamp=2.5e4)
popt = cwresult[0]
pcov = cwresult[1]
perr = np.sqrt(np.diag(pcov))
fit = cwresult[2]
fitg = cwresult[3]
dips = cwresult[4]

plt.close('all')

plt.figure(1,[8,6])

plt.plot(x,y,color='#ED1035',marker='.',label="data")
plt.plot(x,fit,color='#97CC04',linewidth=2.0)
#plt.plot(x,fitg,color='#2D7DD2',linewidth=2.0)
#plt.plot(dips[0][0],dips[1][0],'b.')
#plt.plot(dips[0][1],dips[1][1],'b.')

plt.xlabel(r'$freq \quad (MHz)$')
plt.ylabel(r'$NV \quad PL \quad (counts)$')
plt.tight_layout()
#pylab.savefig('/Users/alec/UCSB/scan_data/images/linecut_lrgunc_'+str(filenum)+'.pdf')

fig = plt.gcf()
fig.canvas.manager.window.raise_()

#print(popt)
bnv = (popt[4]-popt[1])/(2*2.8)
bz = 12.0
print(180*np.arccos(bnv/bz)/np.pi)

