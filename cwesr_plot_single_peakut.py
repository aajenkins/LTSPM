# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""

import numpy as np
import matplotlib.pyplot as plt
import cwesr_fit_single as cwesr
#import math as math
#import matplotlib.gridspec as gridspec
#import matplotlib.gridspec as gridspec

basepath = '/Users/alec/UCSB/scan_data/1566-esrdata/fittrack2'
fitdata = []

filenum = 2852
num_avg = 3

def func(x, *params):
        y = np.zeros_like(x)
        c = params[0]
        for i in range(1, len(params), 3):
            ctr = params[i]
            amp = params[i+1]
            wid = params[i+2]
            #y = y + (amp * (np.exp( -((x - ctr-2.3)/wid)**2)+np.exp( -((x - ctr)/wid)**2)+np.exp( -((x - ctr+2.3)/wid)**2)) + c
            y = y + (-abs(amp * (wid/2)**2)/((x-ctr)**2+(wid/2)**2))
        y=y+c
        return y

plt.close('all')
plt.figure(1,[10,8])
   
filepath = basepath+str(filenum).zfill(6)
data = np.loadtxt(filepath+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]

edata = np.transpose(data)
x, y = edata

cwresult = cwesr.cwesr_fit(x,y)
fit = cwresult[1]
fitg = cwresult[2]
indexes = cwresult[3]

print(cwresult[0][4]-cwresult[0][1])
#plt.plot(x[indexes], y[indexes],'r.')
plt.plot(indexes[0], indexes[1],'r.')
plt.plot(x, y)
plt.plot(x, fit, '-r')
plt.plot(x, fitg, '-b')
#plt.plot(x, backfit, '-g')

#rem = y-fit
#back = [rem[1], 2820, 3e3, 200,3200, 3e3, 200]
#bpopt, bpcov = curve_fit(func, x, rem, p0=back)
#backfit = func(x, *bpopt)
#plt.plot(x,rem)
#plt.plot(x, backfit, '-g')
#print(bpopt)
#fig=plt.gcf()
#fig.canvas.manager.window.activateWindow()
#fig.canvas.manager.window.raise_()