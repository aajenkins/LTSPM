# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:00:36 2017

@author: alec
"""

from scipy.optimize import curve_fit
import peakutils
from peakutils.plot import plot as pplot
import numpy as np
import matplotlib.pyplot as plt
from pylab import get_current_fig_manager

def fit_gaussian(x, *params):
        y = np.zeros_like(x)
        c = params[0]
        for i in range(1, len(params), 3):
            ctr = params[i]
            amp = params[i+1]
            wid = params[i+2]
            #y = y + (amp * (np.exp( -((x - ctr-2.3)/wid)**2)+np.exp( -((x - ctr)/wid)**2)+np.exp( -((x - ctr+2.3)/wid)**2)) + c
            y = y + (-abs(amp * (1/np.sqrt(2*np.pi*(wid**2))) * np.exp(-((x-ctr)**2)/(2*(wid**2)))))
        y=y+c
        return y

dgamp = 1.5e4
dgwidth = 4
lbounds2 = [0,2500,1e3,4,2500,1e3,4]
ubounds2 = [3e5,3200,2e5,50,3200,2e5,50]
pdheight = 1500
maxcenshift = 20
defaultf1 = 2855
defaultf2 = 2886

#import math as math
#import matplotlib.gridspec as gridspec
#import matplotlib.gridspec as gridspec

basepath = '/Users/alec/UCSB/scan_data/1847-esrdata/ff2'
fitdata = []

filenum = 1
num_avg = 2


plt.close('all')
plt.figure(1,[10,8])
   
filepath = basepath+str(filenum).zfill(6)
data = np.loadtxt(filepath+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]

edata = np.transpose(data)
x, y = edata

y = np.add(np.max(y),-y)

indexes = peakutils.indexes(y, thres=0.5, min_dist=3)
print(indexes)
print(x[indexes], y[indexes])


plt.plot(x, y)
plt.plot(x[indexes], y[indexes],'r.')
fig = get_current_fig_manager()
fig.canvas.manager.window.raise_()