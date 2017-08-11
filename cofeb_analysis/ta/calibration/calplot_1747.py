# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:17:32 2016

@author: alec
"""

import time
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

pi = np.pi
#cal parameters
#    theta = 60.4*np.pi/180
bz0 = 11.0
phi = 0*np.pi/180

mstnm = 654577.0
theta = 1.05602865206




#pi = np.pi
#bz0 = 11.0
#phi = 0*pi/180

#file constants
xres = 100
yres = 10
scan_size = 0.8*5000
dres = scan_size/xres

hlist = np.zeros(2*yres)
thetalist = np.zeros(2*yres)
mstlist = np.zeros(2*yres)

ffdata = ls.load_ff('/Users/alec/UCSB/scan_data/1747-esrdata/fitdata.txt',xres,yres,20)

x = np.arange(0,dres*xres,dres)

plt.close('all')

ploth = 1
plt.figure(1,[17,10])
gs = gridspec.GridSpec(5, 4)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)



for j in range(0,10):
    y = ffdata[0][j,:]
    ye = ffdata[1][j,:]
    ry = np.flipud(ffdata[2][j,:])
    rye = np.flipud(ffdata[3][j,:])

    guess = [80, 1600]
    rguess = [80, 1400]
    try:
        popt, pcov = curve_fit(cal_func, x, y, p0=guess)
    except:
        popt = np.zeros(2)
        pcov = np.zeros((2,2))
        print('fit fail')
    try:
        rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess)
    except:
        rpopt = np.zeros(2)
        rpcov = np.zeros((2,2))
        print('fit fail')

    # mstlist[2*j] = popt[0]
    # mstlist[2*j+1] = rpopt[0]
    hlist[2*j] = popt[0]
    hlist[2*j+1] = rpopt[0]
    # thetalist[2*j] = popt[3]
    # thetalist[2*j+1] = rpopt[3]

    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2)])
    plt.plot(x,cal_func(x,*guess),'g-')
    plt.errorbar(x,y,yerr=ye,color='#000000',fmt='.')
    plt.plot(x,cal_func(x,*popt),'r-')
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2+1)])
    plt.plot(x,cal_func(x,*rguess),'g-')
    plt.errorbar(x,ry,yerr=rye,color='#000000',fmt='.')
    plt.plot(x,cal_func(x,*rpopt),'r-')
plt.show()

fig = plt.gcf()
fig.canvas.manager.window.raise_()


print('Ms*t mean = '+str(np.mean(mstlist))+' +/- '+str(np.std(mstlist)))
print('h mean = '+str(np.mean(hlist))+' +/- '+str(np.std(hlist)))
print('theta mean = '+str(np.mean(thetalist))+' +/- '+str(np.std(thetalist)))
