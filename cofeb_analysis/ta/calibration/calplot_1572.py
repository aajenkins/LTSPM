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
import json
import load_scan as lscan
import plotting.format_plots_tkagg as fp

pi = np.pi
#cal parameters
#    theta = 60.4*np.pi/180
bz0 = 70

phi = pi

mstnm = (2e-3)*(1e9)

def cal_func(x, *args):
    bnv = np.zeros_like(x)

    h = args[0]
    x0 = args[1]
    theta = 0.962234639607
    # theta = args[2]
    mstnm = args[3]
    # phi = args[3]
    bx = -(2e-3)*mstnm*(h/(h**2+(x-x0)**2))
    bz = (2e-3)*mstnm*((x-x0)/(h**2+(x-x0)**2))
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+(bz+bz0)*np.cos(theta))
    return bnv

#file constants
xres = 150
yres = 5
scan_size = 1.2*5000
dres = scan_size/xres

mstlist = np.zeros(2*yres)
hlist = np.zeros(2*yres)
thetalist = np.zeros(2*yres)
# philist = np.zeros(2*yres)
r2list = np.zeros(2*yres)

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/1572-esrdata/fitdata.txt',xres,yres,printNVCalcError=True)
x = np.arange(0,dres*xres,dres)

BzMeanEnd = 0
BzStdEnd = 0
BzTails = [0,1,2,3,4]
BzEnd = np.zeros((len(BzTails),12))

for i in range(len(BzTails)):
    BzEnd[i] = np.concatenate((ffdata[0][BzTails[i],-6:], ffdata[2][BzTails[i],:6]))

BzMeanEnd = np.mean(BzEnd)
BzStdEnd = np.std(BzEnd)

thetaMean = np.arccos(BzMeanEnd/bz0)
thetaError = np.sqrt((1/(bz0**2 - BzMeanEnd**2)) + ( (BzMeanEnd**2) / ((bz0**2)*(bz0**2 - BzMeanEnd**2)) ))

print(np.sqrt((1/(bz0**2 - BzMeanEnd**2))))
print(np.sqrt(( (BzMeanEnd**2) / ((bz0**2)*(bz0**2 - BzMeanEnd**2)) )))
print("BzMeanEnd = "+str(BzMeanEnd))
print("thetaMean = "+str(thetaMean))
print("thetaError = "+str(thetaError))

plt.close('all')

ploth = 1
plt.figure(1,[17,10])
gs = gridspec.GridSpec(5, 4)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)

for j in range(0,5):
    y = ffdata[0][j,:]
    ye = ffdata[1][j,:]
    ry = np.flipud(ffdata[2][j,:])
    rye = np.flipud(ffdata[3][j,:])

    ymax = np.max(y)
    yargmax = np.argmax(y)
    ryargmax = np.argmax(ry)
    # yepeak = ye.copy()
    # ryepeak = rye.copy()
    # yepeak[yargmax-2:yargmax+2] = 0.1
    # ryepeak[ryargmax-2:ryargmax+2] = 0.1

    fitLength = 20
    offset = 0
    xShort = np.linspace(-(fitLength-offset)*xres/2, (fitLength+offset)*xres/2, 2*fitLength+1)
    xShortFit = np.linspace(-(fitLength-offset)*xres/2, (fitLength+offset)*xres/2, 10*fitLength+1)
    yShort = y[yargmax-(fitLength-offset):yargmax+(fitLength+offset)+1]
    yeShort = ye[yargmax-(fitLength-offset):yargmax+(fitLength+offset)+1]
    ryShort = ry[ryargmax-(fitLength-offset):ryargmax+(fitLength+offset)+1]
    ryeShort = rye[ryargmax-(fitLength-offset):ryargmax+(fitLength+offset)+1]

    yargmax = 0
    ryargmax = 0
    # yeShort[fitLength-1:fitLength+2] = 0.2
    # ryeShort[fitLength-1:fitLength+2] = 0.2

    guess = [150, yargmax*dres, 1, 0]
    rguess = [150, ryargmax*dres, 1, 0]
    try:
        # popt, pcov = curve_fit(cal_func, x, y, sigma=ye, p0=guess)
        popt, pcov = curve_fit(cal_func, xShort, yShort, sigma=yeShort, p0=guess)
    except:
        popt = np.zeros(4)
        pcov = np.zeros((4,4))
        print('fit fail')
    try:
        # rpopt, rpcov = curve_fit(cal_func, x, ry, sigma=rye, p0=rguess)
        rpopt, rpcov = curve_fit(cal_func, xShort, ryShort, sigma=ryeShort, p0=rguess)
#        rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess, sigma=ryms)
    except:
        rpopt = np.zeros(4)
        rpcov = np.zeros((4,4))
        print('fit fail')
    hlist[2*j] = popt[0]
    hlist[2*j+1] = rpopt[0]
    # philist[2*j] = popt[3]
    # philist[2*j+1] = rpopt[3]
    thetalist[2*j] = popt[2]
    thetalist[2*j+1] = rpopt[2]
    mstlist[2*j] = popt[3]
    mstlist[2*j+1] = rpopt[3]
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

plt.show()

# print('h mean = '+str(np.mean(np.r_[hlist[:5],hlist[[7,9]]]))+' +/- '+str(np.std(np.r_[hlist[:5],hlist[[7,9]]])))
print('h mean = '+str(np.mean(hlist))+' +/- '+str(np.std(hlist)))
print('mst mean = '+str(np.mean(mstlist))+' +/- '+str(np.std(mstlist)))
print('theta mean = '+str(np.mean(thetalist))+' +/- '+str(np.std(thetalist)))
# print('phi mean = '+str(np.mean(philist))+' +/- '+str(np.std(philist)))

#print('bz0 mean = '+str(np.mean(bz0list))+' +/- '+str(np.std(bz0list)))
print("r2 mean = "+str(np.mean(r2list)))
