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
import format_plot as fp

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}

pi = np.pi
#cal parameters
#    theta = 60.4*np.pi/180
bz0 = 12

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
filespec = 'Msfixed'
cal_params = np.loadtxt(path+'cal_parameters_'+filespec+'.txt', delimiter=',')

Ms = cal_params[0]
t = cal_params[1]
phi = cal_params[3]

mstnm = Ms*t*(1e10)

matplotlib.rc('font', **font)

def cal_func(x, *args):
    bnv = np.zeros_like(x)

    h = args[0]
    x0 = args[1]
    theta = 45*pi/180#1.0335227055519
    # theta = args[2] * pi / 180
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
# philist = np.zeros(2*yres)

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/1773-esrdata/fitdata.txt',xres,yres,15)
x = np.arange(0,dres*xres,dres)

plt.close('all')

ploth = 1
plt.figure(1,[17,10])
gs = gridspec.GridSpec(5, 4)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)

BzMeanEnd = 0
BzTails = [0,1,2,3,4,6,7]

for i in BzTails:
    BzMeanEnd += np.mean(ffdata[0][i,-6:])
    BzMeanEnd += np.mean(ffdata[2][i,:6])
BzMeanEnd = BzMeanEnd/(2*len(BzTails))

for j in range(0,10):
    y = ffdata[0][j,:]
    ye = ffdata[1][j,:]
    ry = np.flipud(ffdata[2][j,:])
    rye = np.flipud(ffdata[3][j,:])

    ymax = np.max(y)
    yargmax = np.argmax(y)
    ryargmax = np.argmax(ry)
    yepeak = ye.copy()
    ryepeak = rye.copy()
    yepeak[yargmax-1:yargmax+1] = 0.1
    ryepeak[ryargmax-1:ryargmax+1] = 0.1

    fitLength = 4
    offset = 1
    xShort = np.linspace(-(fitLength-offset)*xres/2, (fitLength+offset)*xres/2, 2*fitLength+1)
    xShortFit = np.linspace(-(fitLength-offset)*xres/2, (fitLength+offset)*xres/2, 10*fitLength+1)
    yShort = y[yargmax-(fitLength-offset):yargmax+(fitLength+offset)+1]
    yeShort = ye[yargmax-(fitLength-offset):yargmax+(fitLength+offset)+1]
    ryShort = ry[ryargmax-(fitLength-offset):ryargmax+(fitLength+offset)+1]
    ryeShort = rye[ryargmax-(fitLength-offset):ryargmax+(fitLength+offset)+1]
    # yeShort[fitLength-1:fitLength+2] = 0.2
    # ryeShort[fitLength-1:fitLength+2] = 0.2

    guess = [50, 0]
    rguess = [50, 0]
    try:
        popt, pcov = curve_fit(cal_func, xShort, yShort, sigma=yeShort, p0=guess)
#        popt, pcov = curve_fit(cal_func, x, y, p0=guess, sigma=yms)
    except:
        popt = np.zeros(4)
        pcov = np.zeros((4,4))
        print('fit fail')
    try:
        rpopt, rpcov = curve_fit(cal_func, xShort, ryShort, sigma=ryeShort, p0=rguess)
#        rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess, sigma=ryms)
    except:
        rpopt = np.zeros(4)
        rpcov = np.zeros((4,4))
        print('fit fail')
    hlist[2*j] = popt[0]
    hlist[2*j+1] = rpopt[0]
    # philist[2*j] = popt[2]
    # philist[2*j+1] = rpopt[2]
    # thetalist[2*j] = popt[2]
    # thetalist[2*j+1] = rpopt[2]
#    bz0list[2*j] = popt[4]
#    bz0list[2*j+1] = rpopt[4]

    plt.figure(1)
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2)])
    plt.plot(xShortFit,cal_func(xShortFit,*guess),'g-')
    plt.errorbar(xShort,yShort,yerr=yeShort,color='#000000',fmt='.')
    plt.plot(xShortFit,cal_func(xShortFit,*popt),'r-')
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2+1)])
    plt.plot(xShortFit,cal_func(xShortFit,*rguess),'g-')
    plt.errorbar(xShort,ryShort,yerr=ryeShort,color='#000000',fmt='.')
    plt.plot(xShortFit,cal_func(xShortFit,*rpopt),'r-')

fp.format_plot(plt, 1200, 900, 0, 50, tight=False)
plt.show()

print('h mean = '+str(np.mean(hlist))+' +/- '+str(np.std(hlist)))
# print('theta mean = '+str(np.mean(thetalist))+' +/- '+str(np.std(thetalist)))
# print('phi mean = '+str(np.mean(philist))+' +/- '+str(np.std(philist)))

#print('bz0 mean = '+str(np.mean(bz0list))+' +/- '+str(np.std(bz0list)))
