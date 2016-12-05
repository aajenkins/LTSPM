# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:17:32 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import numpy as np
import math as math
import load_scan as ls

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

def cal_func(x, *params):
        bnv = np.zeros_like(x)
        
        mst = params[0]
        h = params[1]
        x0 = params[2]
        theta = params[3]*pi/180
#        phi = params[4]*pi/180
        
        bx = (2e-3)*mst*(h/(h**2+(x-x0)**2))
        bz = (2e-3)*mst*((x-x0)/(h**2+(x-x0)**2))
        #y = y + (amp * (np.exp( -((x - ctr-2.3)/wid)**2)+np.exp( -((x - ctr)/wid)**2)+np.exp( -((x - ctr+2.3)/wid)**2)) + c
        bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+(bz+bz0)*np.cos(theta))
        return bnv
        
#material parameters
pi = np.pi
sres = 10
bz0 = 11.0
phi = 0

#file constants
dres = 5000*0.8/150

xres = 120
yres = 5

ffdata = ls.load_ff('/Users/alec/UCSB/scan_data/1747-esrdata/fitdata.txt',xres,yres,7)

x = np.arange(0,dres*150,dres)

plt.close('all')
    
ploth = 1
plt.figure(1,[17,10])
gs = gridspec.GridSpec(5, 4)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)
    
for j in range (0,20):
    y = ffdata[0][j,:]
    ye = ffdata[1][j,:]
    if np.mod(j,2) == 1:
        y = np.flipud(y)
        ye = np.flipud(ye)
    
    guess = [2e6, 130, 3200, 54]
    try:
        popt, pcov = curve_fit(cal_func, x, y, p0=guess)
    except:
        popt = np.zeros(4)
        pcov = np.zeros((4,4))   
    
    csubplot = plt.subplot(gs[(j%5),math.floor(j/5)])
    #plt.plot(x,cal_func(x,*guess),'r-')
    plt.errorbar(x,y,yerr=ye,color='#000000',fmt='.')
    plt.plot(x,cal_func(x,*popt),'r-')
    
plt.show()

fig = plt.gcf()
fig.canvas.manager.window.raise_()
