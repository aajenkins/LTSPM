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
import load_scan as lscan

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}

pi = np.pi
#cal parameters
#    theta = 60.4*np.pi/180
bz0 = 12
phi = 0*np.pi/180
mst = 6.22e-4
mstnm = mst/(1.0e-9)

matplotlib.rc('font', **font)

def cal_func(x, *args):
    bnv = np.zeros_like(x)

    h = args[0]
    x0 = args[1]
    theta = args[2] * pi / 180
    bx = (2e-3)*mstnm*(h/(h**2+(x-x0)**2))
    bz = -(2e-3)*mstnm*((x-x0)/(h**2+(x-x0)**2))
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+(bz+bz0)*np.cos(theta))
    return bnv

#file constants
xres = 100
yres = 10
scan_size = 0.8*5000
dres = scan_size/xres

hlist = np.zeros(2*yres)
thetalist = np.zeros(2*yres)

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/1773-esrdata/fitdata.txt',xres,yres,15)
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

    guess = [70, 2200, 55]
    rguess = [70, 2000, 55]
    try:
        popt, pcov = curve_fit(cal_func, x, y, sigma=ye, p0=guess)
#        popt, pcov = curve_fit(cal_func, x, y, p0=guess, sigma=yms)
    except:
        popt = np.zeros(4)
        pcov = np.zeros((4,4))
        print('fit fail')
    try:
        rpopt, rpcov = curve_fit(cal_func, x, ry, sigma=rye, p0=rguess)
#        rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess, sigma=ryms)
    except:
        rpopt = np.zeros(4)
        rpcov = np.zeros((4,4))
        print('fit fail')
    hlist[2*j] = popt[0]
    hlist[2*j+1] = rpopt[0]
    thetalist[2*j] = popt[2]
    thetalist[2*j+1] = rpopt[2]
#    bz0list[2*j] = popt[4]
#    bz0list[2*j+1] = rpopt[4]
#

    plt.figure(1)
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2)])
    plt.plot(x,cal_func(x,*guess),'g-')
    plt.errorbar(x,y,yerr=ye,color='#000000',fmt='.')
    plt.plot(x,cal_func(x,*popt),'r-')
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

print('h mean = '+str(np.mean(hlist))+' +/- '+str(np.std(hlist)))
print('theta mean = '+str(np.mean(thetalist))+' +/- '+str(np.std(thetalist)))
#print('bz0 mean = '+str(np.mean(bz0list))+' +/- '+str(np.std(bz0list)))
