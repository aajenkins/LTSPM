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
from scipy.integrate import quad
import numpy as np
import json
import load_scan as lscan
import plotting.format_plots_tkagg as fp


pi = np.pi


# def height_calibration(scannum):

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
DMI = material_params['DMI']
DWWidth = material_params['DWWidth']
Aex = material_params['Aex']

theta = scan_params['theta']
thetaError = scan_params['thetaError']
phi = scan_params['phi']
calibrationScanNum = scan_params['calibrationScanNum']
xresCal = scan_params['xresCal']
yresCal = scan_params['yresCal']
scanSizeCal = scan_params['scanSizeCal']
zfieldCal = scan_params['zfieldCal']*(1e4)

dres = scanSizeCal/xresCal
xi = 2*Aex/DMI
mu0 = 4*pi*(1e-7)
B0 = (1e4)*mu0*Ms*t/(2*pi)
xmcut = 5e-6#100*DWWidth

def bx_tilt(xm, x, h, x0):
    mx = (DWWidth/xi) * np.exp((xm-x0)/DWWidth)
    mz = np.sqrt(1-mx**2)
    return -( mx/((x-xm)**2 + h**2) - 2*(x-xm)*(mx*(x-xm) + mz*h)/(((x-xm)**2 + h**2)**2) )


def bz_tilt(xm, x, h, x0):
    mx = (DWWidth/xi) * np.exp((xm-x0)/DWWidth)
    mz = np.sqrt(1-mx**2)
    return -( mz/((x-xm)**2 + h**2) - 2*h*(mx*(x-xm) + mz*h)/(((x-xm)**2 + h**2)**2) )

def tilt_edge(x, *params):
    h = params[0]
    x0 = params[1]

    Bx0 = B0 * h/(((x+xmcut)**2) + h**2)
    Bz0 = -B0 * (x+xmcut)/(((x+xmcut)**2) + h**2)

    Bx = Bx0 + B0 * (quad(bx_tilt, x0-xmcut, x0, args=(x,h,x0))[0])
    Bz = Bz0 + B0 * (quad(bz_tilt, x0-xmcut, x0, args=(x,h,x0))[0])

    BNV = np.abs(Bx*np.sin(theta)*np.cos(phi)+(Bz+zfieldCal)*np.cos(theta))

    return BNV

vtilt_edge = np.vectorize(tilt_edge)

def edge_line(x, *params):
    bnv = np.zeros_like(x)

    h = params[0]
    x0 = params[1]

    Bliner = B0 * (DWWidth/xi) / (np.sqrt((x-x0)**2 + h**2))
    Blinex = Bliner*((x-x0)/np.sqrt((x-x0)**2 + h**2))
    Blinez = Bliner*(h/np.sqrt((x-x0)**2 + h**2))

    Bx = Blinex + B0*(h/(h**2+(x-x0)**2))
    Bz = Blinez - B0*((x-x0)/(h**2+(x-x0)**2))
    BNV = np.abs(Bx*np.sin(theta)*np.cos(phi)+(Bz+zfieldCal)*np.cos(theta))

    return BNV

def edge(x, *params):
    bnv = np.zeros_like(x)

    h = params[0]
    x0 = params[1]

    bx = B0*(h/(h**2+(x-x0)**2))
    bz = -B0*((x-x0)/(h**2+(x-x0)**2))
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+(bz+zfieldCal)*np.cos(theta))

    return bnv


cal_func = vtilt_edge

hlist = np.zeros(2*yresCal)
# mstlist = np.zeros(2*yresCal)
r2list = np.zeros(2*yresCal)

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(calibrationScanNum)+'-esrdata/fitdata.txt',
xresCal,yresCal,maxfgrad=20, printNVCalcError=True)
x = np.arange(0,dres*xresCal,dres)

# xt = np.arange(-2*dres*xresCal,2*dres*xresCal,dres)
# bnv = cal_func(xt,6e-8,0)

# fig, ax = plt.subplots()
# plt.plot(xt*(1e6), bnv)
# plt.show()

BzMeanEnd = 0
BzStdEnd = 0
BzTails = [0,1,2,3,4,6,7]
BzEnd = np.zeros((len(BzTails),12))

for i in range(len(BzTails)):
    BzEnd[i] = np.concatenate((ffdata[0][BzTails[i],-6:], ffdata[2][BzTails[i],:6]))

BzMeanEnd = np.mean(BzEnd)
BzStdEnd = np.std(BzEnd)

thetaMean = np.arccos(BzMeanEnd/zfieldCal)
thetaError = np.sqrt((1/(zfieldCal**2 - BzMeanEnd**2)) + ( (BzMeanEnd**2) / ((zfieldCal**2)*(zfieldCal**2 - BzMeanEnd**2)) ))

print(np.sqrt((1/(zfieldCal**2 - BzMeanEnd**2))))
print(np.sqrt(( (BzMeanEnd**2) / ((zfieldCal**2)*(zfieldCal**2 - BzMeanEnd**2)) )))
print("BzMeanEnd = "+str(BzMeanEnd))
print("thetaMean = "+str(thetaMean))
print("thetaError = "+str(thetaError))

fitLength = 20
offset = 0
xShort = np.linspace(-(fitLength-offset)*dres, (fitLength+offset)*dres, 2*fitLength+1)
xShortFit = np.linspace(-(fitLength-offset)*dres, (fitLength+offset)*dres, 10*fitLength+1)
xFit = np.linspace(np.min(xShortFit), np.max(xShortFit), 5*len(x))

fits = np.empty((yresCal), dtype=object)
rfits = np.empty((yresCal), dtype=object)
yShorts = np.empty((yresCal), dtype=object)
yeShorts = np.empty((yresCal), dtype=object)
ryShorts = np.empty((yresCal), dtype=object)
ryeShorts = np.empty((yresCal), dtype=object)

plt.close('all')

fig, ax = plt.subplots(figsize=(12,7))
gs = gridspec.GridSpec(5, 4)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)

for j in range(0,yresCal):

    print(j)

    y = ffdata[0][j,:]
    ye = ffdata[1][j,:]
    ry = np.flipud(ffdata[2][j,:])
    rye = np.flipud(ffdata[3][j,:])

    yargmax = np.argmax(y)
    ryargmax = np.argmax(ry)

    yShorts[j] = y[yargmax-(fitLength-offset):yargmax+(fitLength+offset)+1]
    yeShorts[j] = ye[yargmax-(fitLength-offset):yargmax+(fitLength+offset)+1]
    ryShorts[j] = ry[ryargmax-(fitLength-offset):ryargmax+(fitLength+offset)+1]
    ryeShorts[j] = rye[ryargmax-(fitLength-offset):ryargmax+(fitLength+offset)+1]

    x = xShort
    y = yShorts[j]
    ye = yeShorts[j]
    ry = ryShorts[j]
    rye = ryeShorts[j]

    guess = [6e-8, 0]
    rguess = [6e-8, 0]
    try:
        # popt, pcov = curve_fit(vtilt_edge, x, y, p0=guess)
        popt, pcov = curve_fit(cal_func, x, y, p0=guess)
    except:
        popt = np.zeros(2)
        pcov = np.zeros((2,2))
        print('fit fail')
    try:
        # rpopt, rpcov = curve_fit(vtilt_edge, x, ry, p0=rguess)
        rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess)
    except:
        rpopt = np.zeros(4)
        rpcov = np.zeros((4,4))
        print('fit fail reverse')
    hlist[2*j] = popt[0]
    hlist[2*j+1] = rpopt[0]
    fits[j] = popt
    rfits[j] = rpopt
    # mstlist[2*j] = popt[2]
    # mstlist[2*j+1] = rpopt[2]

    ss_res = np.sum((y - cal_func(x,*popt)) ** 2)
    ss_tot = np.sum((y - np.mean(y) ** 2))
    r2 = 1 - (ss_res / ss_tot)
    rss_res = np.sum((ry - cal_func(x,*rpopt)) ** 2)
    rss_tot = np.sum((ry - np.mean(ry)) ** 2)
    rr2 = 1 - (ss_res / ss_tot)

    r2list[2*j] = r2
    r2list[2*j+1] = rr2

    plt.figure(1)
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2)])
    plt.plot(xFit,cal_func(xFit,*guess))
    plt.errorbar(x,y,yerr=ye,color='#000000',fmt='.')
    plt.plot(xFit,cal_func(xFit,*popt))
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2+1)])
    plt.plot(xFit,cal_func(xFit,*rguess))
    plt.errorbar(x,ry,yerr=rye,color='#000000',fmt='.')
    plt.plot(xFit,cal_func(xFit,*rpopt))

fp.format_plots(plt, small=False, tight=False)

snum = 7
xFitum = xFit*(1e6)
xum=x*(1e6)
fig, ax = plt.subplots(figsize=[6,4])
plt.plot(xFitum,cal_func(xFit,*rfits[snum]))
plt.errorbar(xum,ryShorts[snum],yerr=ryeShorts[snum],
             color='#000000',fmt='.')
plt.xlabel(r'x $\mu$m')
plt.ylabel(r'B$_{NV}$ (G)')
plt.subplots_adjust(left=0.1, right=1.0, bottom=0.13, top=1.0)

plt.show()


print('h mean = '+str(np.mean(hlist))+' +/- '+str(np.std(hlist)))
# print('mst mean = '+str(np.mean(mstlist))+' +/- '+str(np.std(mstlist)))
print("r2 mean = "+str(np.mean(r2list)))

# if __name__ == "__main__":
#     import sys
#     if (len(sys.argv) == 2):
#         height_calibration(int(sys.argv[1]))
#     else:
#         print('enter scan number')
