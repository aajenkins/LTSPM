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
import load_scan as ls

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}

pi = np.pi
#cal parameters
#    theta = 60.4*np.pi/180
bz0 = 12
phi = 0*np.pi/180
    
matplotlib.rc('font', **font)

def cal_func(x, *params):
    bnv = np.zeros_like(x)
        
    mst = params[0]
    h = params[1]
    x0 = params[2]
    theta = params[3]*pi/180
#    theta = 61.507*pi/180
#    bz0 = params[4]
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
mstlist = np.zeros(2*yres)
#bz0list = np.zeros(2*yres)

ffdata = ls.load_ff('/Users/alec/UCSB/scan_data/1773-esrdata/fitdata.txt',xres,yres,10)
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
    
    ymax = np.max(y)
    yargmax = np.argmax(y)
    ryargmax = np.argmax(ry)
    yms = np.zeros(xres)
    yms.fill(5)
    yms[yargmax:yargmax+10] = 1
    ryms = np.zeros(xres)
    ryms.fill(5)
    ryms[ryargmax:ryargmax+10] = 1
    
    guess = [8e5, 80, 2200, 55]
    rguess = [8e5, 80, 2000, 55]
    try:
        popt, pcov = curve_fit(cal_func, x, y, p0=guess)
#        popt, pcov = curve_fit(cal_func, x, y, p0=guess, sigma=yms)
    except:
        popt = np.zeros(4)
        pcov = np.zeros((4,4))
        print('fit fail')
    try:
        rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess)
#        rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess, sigma=ryms)
    except:
        rpopt = np.zeros(4)
        rpcov = np.zeros((4,4)) 
        print('fit fail')
    mstlist[2*j] = popt[0]
    mstlist[2*j+1] = rpopt[0]
    hlist[2*j] = popt[1]
    hlist[2*j+1] = rpopt[1]
    thetalist[2*j] = popt[3]
    thetalist[2*j+1] = rpopt[3]
#    bz0list[2*j] = popt[4]
#    bz0list[2*j+1] = rpopt[4]
#
        
    plt.figure(1)
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2)])
    plt.plot(x,cal_func(x,*guess),'g-')
    plt.errorbar(x,y,yerr=ye,color='#000000',fmt='.')
    plt.plot(x,cal_func(x,*popt),'r-')
    plt.plot(yargmax*dres,ymax,'co')
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2+1)])
    plt.plot(x,cal_func(x,*rguess),'g-')
    plt.errorbar(x,ry,yerr=rye,color='#000000',fmt='.')
    plt.plot(x,cal_func(x,*rpopt),'r-')
plt.show()
#


fig = plt.gcf()
fig.canvas.manager.window.raise_()

spnum = 0
dstart = 10
dend = 90
y = ffdata[0][spnum,dstart:dend]
ye = ffdata[1][spnum,dstart:dend]
ry = np.flipud(ffdata[2][spnum,:])[dstart:dend]
rye = np.flipud(ffdata[3][spnum,:])[dstart:dend]
xc = x[dstart:dend]
lc = len(xc)

ymax = np.max(y)
yargmax = np.argmax(y)
ryargmax = np.argmax(ry)
yms = np.zeros(lc)
yms.fill(5)
yms[yargmax:yargmax+1] = 1
ryms = np.zeros(lc)
ryms.fill(5)
ryms[ryargmax:ryargmax+10] = 1
    
try:
#     popt, pcov = curve_fit(cal_func, x, y, p0=guess)
     popt, pcov = curve_fit(cal_func, xc, y, p0=guess, sigma=yms)
#    rpopt, rpcov = curve_fit(cal_func, xc, ry, p0=rguess)
#    rpopt, rpcov = curve_fit(cal_func, xc, ry, p0=rguess, sigma=ryms)
except:
    popt = np.zeros(4)
    pcov = np.zeros((4,4)) 
    print('fit fail')

plt.figure(2,[8,6])
plt.errorbar(xc,y,yerr=ye,color='#ED1035', fmt='.')
plt.plot(xc,cal_func(xc,*popt),color='#2D7DD2',linewidth=2.0)
plt.xlabel(r'$x \quad (nm)$')
plt.ylabel(r'$B_{NV} \quad (G)$')
plt.tight_layout()
#pylab.savefig('/Users/alec/UCSB/scan_data/images/calibration_ex_'+str(filenum)+'.png')

fig = plt.gcf()
fig.canvas.manager.window.raise_()


    
print('Ms*t mean = '+str(np.mean(mstlist))+' +/- '+str(np.std(mstlist)))
print('h mean = '+str(np.mean(hlist))+' +/- '+str(np.std(hlist)))
print('theta mean = '+str(np.mean(thetalist))+' +/- '+str(np.std(thetalist)))
#print('bz0 mean = '+str(np.mean(bz0list))+' +/- '+str(np.std(bz0list)))

