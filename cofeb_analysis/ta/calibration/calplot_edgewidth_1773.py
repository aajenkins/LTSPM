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

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}

pi = np.pi
#cal parameters
#    theta = 60.4*np.pi/180

scannum = 1760

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
scan_path = path+str(scannum)+'/'
scan_params_path = scan_path+'scan_parameters.json'
material_params_path = path+'material_parameters.json'
with open(scan_params_path, 'r') as fread:
    scan_params = json.load(fread)
with open(material_params_path, 'r') as fread:
    material_params = json.load(fread)

Ms = material_params['Ms']
t = material_params['t']
Msterror = material_params['MstError']

theta = scan_params['theta']
thetaError = scan_params['thetaError']
phi = scan_params['phi']
calibrationScanNum = scan_params['calibrationScanNum']
xresCal = scan_params['xresCal']
yresCal = scan_params['yresCal']
scanSizeCal = scan_params['scanSizeCal']*(1e9)
zfieldCal = scan_params['zfieldCal']*(1e4)

mstnm = Ms*t*(1e9)
dres = scanSizeCal/xresCal

dwWidth = 3e-8

matplotlib.rc('font', **font)

def cal_func(x, *args):
    bnv = np.zeros_like(x)

    h = args[0]
    x0 = args[1]
    dwWidthFit = args[2]
    # theta = args[2] * pi / 180
    N = 20
    jlist = np.arange(N)
    bdenomlist = np.zeros_like(x)
    for i in range(len(x)):
        bdenomlist[i] = (1/N)*np.sum((1/(h**2+(x[i]+dwWidthFit*np.arctanh(jlist/N)-x0)**2)))
    bx = (2e-3)*mstnm*h*bdenomlist
    bz = -(2e-3)*mstnm*(x-x0)*bdenomlist
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+(bz+zfieldCal)*np.cos(theta))
    return bnv


hlist = np.zeros(2*yresCal)
thetalist = np.zeros(2*yresCal)
# philist = np.zeros(2*yresCal)
dwlist = np.zeros(2*yresCal)
r2list = np.zeros(2*yresCal)

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(calibrationScanNum)+'-esrdata/fitdata.txt',
xresCal,yresCal,maxfgrad=20, printNVCalcError=True)
x = np.arange(0,dres*xresCal,dres)

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

fits = np.empty((yresCal), dtype=object)
rfits = np.empty((yresCal), dtype=object)
yShorts = np.empty((yresCal), dtype=object)
yeShorts = np.empty((yresCal), dtype=object)
ryShorts = np.empty((yresCal), dtype=object)
ryeShorts = np.empty((yresCal), dtype=object)

fitLength = 25
offset = 0
xShort = np.linspace(-(fitLength-offset)*xresCal/2, (fitLength+offset)*xresCal/2, 2*fitLength+1)
xShortFit = np.linspace(-(fitLength-offset)*xresCal/2, (fitLength+offset)*xresCal/2, 10*fitLength+1)
xFit = np.linspace(np.min(xShortFit), np.max(xShortFit), 5*len(x))

for j in range(0,10):
    y = ffdata[0][j,:]
    ye = ffdata[1][j,:]
    ry = np.flipud(ffdata[2][j,:])
    rye = np.flipud(ffdata[3][j,:])

    ymax = np.max(y)
    yargmax = np.argmax(y)+1
    ryargmax = np.argmax(ry)+1

    yShorts[j] = y[yargmax-(fitLength-offset):yargmax+(fitLength+offset)+1]
    yeShorts[j] = ye[yargmax-(fitLength-offset):yargmax+(fitLength+offset)+1]
    ryShorts[j] = ry[ryargmax-(fitLength-offset):ryargmax+(fitLength+offset)+1]
    ryeShorts[j] = rye[ryargmax-(fitLength-offset):ryargmax+(fitLength+offset)+1]

    yargmax = 0
    ryargmax = 0

    x = xShort
    y = yShorts[j]
    ye = yeShorts[j]
    ry = ryShorts[j]
    rye = ryeShorts[j]

    guess = [60, yargmax*dres, dwWidth]
    rguess = [60, ryargmax*dres, dwWidth]
    try:
        # popt, pcov = curve_fit(cal_func, x, y, sigma=ye, p0=guess)
        popt, pcov = curve_fit(cal_func, x, y, sigma=ye, p0=guess)
#        popt, pcov = curve_fit(cal_func, x, y, p0=guess, sigma=yms)
    except:
        popt = np.zeros(4)
        pcov = np.zeros((4,4))
        print('fit fail')
    try:
        rpopt, rpcov = curve_fit(cal_func, x, ry, sigma=rye, p0=rguess)
        # rpopt, rpcov = curve_fit(cal_func, x, ry, sigma=rye, p0=rguess)
#        rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess, sigma=ryms)
    except:
        rpopt = np.zeros(4)
        rpcov = np.zeros((4,4))
        print('fit fail')
    hlist[2*j] = popt[0]
    hlist[2*j+1] = rpopt[0]
    dwlist[2*j] = popt[2]
    dwlist[2*j+1] = rpopt[2]
    fits[j] = popt
    rfits[j] = rpopt

    ss_res = np.sum((y - cal_func(x,*popt)) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot)
    rss_res = np.sum((ry - cal_func(x,*rpopt)) ** 2)
    rss_tot = np.sum((ry - np.mean(ry)) ** 2)
    rr2 = 1 - (ss_res / ss_tot)

    r2list[2*j] = r2
    r2list[2*j+1] = rr2
    # philist[2*j] = popt[2]
    # philist[2*j+1] = rpopt[2]
    # thetalist[2*j] = popt[2]
    # thetalist[2*j+1] = rpopt[2]
#    bz0list[2*j] = popt[4]
#    bz0list[2*j+1] = rpopt[4]

    plt.figure(1)
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2)])
    plt.plot(xFit,cal_func(xFit,*guess),'g-')
    plt.errorbar(x,y,yerr=ye,color='#000000',fmt='.')
    plt.plot(xFit,cal_func(xFit,*popt),'r-')
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2+1)])
    plt.plot(xFit,cal_func(xFit,*rguess),'g-')
    plt.errorbar(x,ry,yerr=rye,color='#000000',fmt='.')
    plt.plot(xFit,cal_func(xFit,*rpopt),'r-')

snum = 9
fig, ax = plt.subplots()
plt.plot(xFit,cal_func(xFit,*fits[snum]))
plt.errorbar(x,yShorts[snum],yerr=yeShorts[snum],
             color='#000000',fmt='.')

fp.format_plots(plt, small=False)
plt.show()

print('h mean = '+str(np.mean(hlist))+' +/- '+str(np.std(hlist)))
# print('theta mean = '+str(np.mean(thetalist))+' +/- '+str(np.std(thetalist)))
# print('phi mean = '+str(np.mean(philist))+' +/- '+str(np.std(philist)))
print('edge width mean = '+str(np.mean(dwlist))+' +/- '+str(np.std(dwlist)))

print("r2 mean = "+str(np.mean(r2list)))

#print('bz0 mean = '+str(np.mean(bz0list))+' +/- '+str(np.std(bz0list)))
