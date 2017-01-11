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
bz0 = 12.0
phi = 0*np.pi/180
    
matplotlib.rc('font', **font)

def cal_func(x, *params):
    bnv = np.zeros_like(x)
        
    mst = params[0]
    h = params[1]
    x0 = params[2]
    theta = params[3]*pi/180
#        phi = params[4]*pi/180
        
    bx = (2e-3)*mst*(h/(h**2+(x-x0)**2))
    bz = -(2e-3)*mst*((x-x0)/(h**2+(x-x0)**2))
    #y = y + (amp * (np.exp( -((x - ctr-2.3)/wid)**2)+np.exp( -((x - ctr)/wid)**2)+np.exp( -((x - ctr+2.3)/wid)**2)) + c
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+(bz+bz0)*np.cos(theta))
    return bnv
        

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
xlist = np.zeros(2*yres)
mstlist = np.zeros(2*yres)

ffdata = ls.load_ff('/Users/alec/UCSB/scan_data/1773-esrdata/fitdata.txt',xres,yres,10)
x = np.arange(0,dres*xres,dres)
    
plt.close('all')
plt.figure(1,[17,10])
gs = gridspec.GridSpec(5, 4)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)



for j in range(0,10):
    y = ffdata[0][j,:]
    ye = ffdata[1][j,:]
    ry = np.flipud(ffdata[2][j,:])
    rye = np.flipud(ffdata[3][j,:])
    
    guess = [8e5, 65, 2200,55]
    rguess = [8e5, 65, 2000,55]
    try:
        popt, pcov = curve_fit(cal_func, x, y, p0=guess)
    except:
        popt = np.zeros(4)
        pcov = np.zeros((4,4))
        print('fit fail')
    try:
        rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess)
    except:
        rpopt = np.zeros(4)
        rpcov = np.zeros((4,4)) 
        print('fit fail')
        
    mstlist[2*j] = popt[0]
    mstlist[2*j+1] = rpopt[0]
    hlist[2*j] = popt[1]
    hlist[2*j+1] = rpopt[1]
    xlist[2*j] = int(np.round(popt[2]/dres))
    xlist[2*j+1] = int(np.round(rpopt[2]/dres))
    thetalist[2*j] = popt[3]
    thetalist[2*j+1] = rpopt[3]
    
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

print(hlist)
print('Ms*t mean = '+str(np.mean(mstlist))+' +/- '+str(np.std(mstlist)))
print('h mean = '+str(np.mean(hlist))+' +/- '+str(np.std(hlist)))
print('theta mean = '+str(np.mean(thetalist))+' +/- '+str(np.std(thetalist)))
xsmax = np.max(xlist)
xsmin = np.min(xlist)

yc = np.zeros(int(100-(xsmax-xsmin)))
xc = x[0:91]

plt.figure(2,[6,6])

for j in range(0,10):
    istart = int(xlist[2*j]-xsmin)
    iend = int(100-(xsmax-xlist[2*j]))
    y = ffdata[0][j,istart:iend]
    ristart = int(xlist[2*j+1]-xsmin)
    riend = int(100-(xsmax-xlist[2*j+1]))
    ry = np.flipud(ffdata[2][j,:])[ristart:riend]
    yc = np.add(yc,y)
    yc = np.add(yc,ry)
    plt.plot(xc,np.add(y,(2*j)*20),'r.')
    plt.plot(xc,np.add(ry,(2*j+1)*20),'b.')

fig = plt.gcf()
fig.canvas.manager.window.raise_()

yc = np.multiply(yc,1/20)


guess = [8e5, 60, 1900,55]
try:
    popt, pcov = curve_fit(cal_func, xc, yc, p0=guess)
except:
    popt = np.zeros(4)
    pcov = np.zeros((4,4))
    print('fit fail')

print(popt[1])

plt.figure(3,[6,6])
plt.plot(xc,yc,'b.')
plt.plot(xc,cal_func(xc,*popt),'r-')

fig = plt.gcf()
fig.canvas.manager.window.raise_()

