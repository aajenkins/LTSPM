# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:17:32 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import numpy as np
import json
import load_scan as lscan
import format_plot as fp

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}

pi = np.pi
#cal parameters
#    theta = 60.4*np.pi/180
bz0 = 40
phi = 180*np.pi/180
mst = 1.07365e-3
mstnm = mst/(1e-9)
sinthetat = -0.047763346911198365
# theta = 0.985437137418

matplotlib.rc('font', **font)

# with open(cal_params_path, 'r') as fread:
#     cal_params = json.load(fread)
#
# Ms = cal_params['Ms']
# t = cal_params['t']
# Msterror = cal_params['Msterror']
# phi = cal_params['phi']

def cal_func(x, *args):
    bnv = np.zeros_like(x)
    h = args[0]
    x0 = args[1]
    theta = args[2]
    # mstv = args[3]
    # theta = 0.983671933677
    bliner = (2e-3)*mstnm*sinthetat/(2*pi*np.sqrt(x**2 + h**2))
    blinex = bliner*(x/np.sqrt(x**2 + h**2))
    blinez = bliner*(h/np.sqrt(x**2 + h**2))
    bx = (2e-3)*mstnm*(h/(h**2+(x-x0)**2))+blinex
    bz = -(2e-3)*mstnm*((x-x0)/(h**2+(x-x0)**2))+blinez
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+(bz+bz0)*np.cos(theta))
    return bnv


#pi = np.pi
#bz0 = 11.0
#phi = 0*pi/180

#file constants
filenum = 1914
xres = 150
yres = 20
scan_size = 0.6*5000
dres = scan_size/xres

hlist = np.zeros(2*yres)
thetalist = np.zeros(2*yres)
# mstlist = np.zeros(2*yres)
#bz0list = np.zeros(2*yres)

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',xres,yres,maxfgrad=25)
x = np.arange(0,dres*xres,dres)

BzMeanEnd = 0
BzStdEnd = 0
BzTails = np.arange(0,yres)
BzEnd = np.zeros((len(BzTails),10))

for i in range(len(BzTails)):
    BzEnd[i] = np.concatenate((ffdata[0][BzTails[i],-5:], ffdata[2][BzTails[i],:5]))

BzMeanEnd = np.mean(BzEnd)
BzStdEnd = np.std(BzEnd)

thetaMean = np.arccos(BzMeanEnd/bz0)
thetaError = np.sqrt((1/(bz0**2 - BzMeanEnd**2)) + ( (BzMeanEnd**2) / ((bz0**2)*(bz0**2 - BzMeanEnd**2)) ))
print(np.sqrt((1/(bz0**2 - BzMeanEnd**2))))
print(np.sqrt(( (BzMeanEnd**2) / ((bz0**2)*(bz0**2 - BzMeanEnd**2)) )))
print('thetaMean = '+str(thetaMean))
print('thetaError = '+str(thetaError))

ploth = 1
fig1, ax1 = plt.subplots()
gs = gridspec.GridSpec(5, 8)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)

guess = [120, 1500, 0.9836]
rguess = guess
rguess[1] = guess[1]-100

for j in range(0,yres):
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

    try:
        popt, pcov = curve_fit(cal_func, x, y, sigma=ye,p0=guess)
    #       popt, pcov = curve_fit(cal_func, x, y, p0=guess, sigma=yms)
    except:
        popt = np.zeros(4)
        pcov = np.zeros((4,4))
        print('fit fail')
    try:
        rpopt, rpcov = curve_fit(cal_func, x, ry, sigma=rye, p0=rguess)
    #       rpopt, rpcov = curve_fit(cal_func, x, ry, p0=rguess, sigma=ryms)
    except:
        rpopt = np.zeros(4)
        rpcov = np.zeros((4,4))
        print('fit fail')

    thetalist[2*j] = popt[2]
    thetalist[2*j+1] = rpopt[2]
    # mstlist[2*j] = popt[3]
    # mstlist[2*j+1] = rpopt[3]
	#    bz0list[2*j] = popt[4]
	#    bz0list[2*j+1] = rpopt[4]
	#
    hlist[2*j] = popt[0]
    hlist[2*j+1] = rpopt[0]
    plt.figure(fig1.number)
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2)])
    plt.plot(x,cal_func(x,*guess),'g-')
    plt.errorbar(x,y,yerr=ye,color='#000000',fmt='.')
    plt.plot(x,cal_func(x,*popt),'r-')
    # plt.plot(yargmax*dres,ymax,'co')
    csubplot = plt.subplot(gs[(j%5),int(np.floor(j/5)*2+1)])
    # plt.plot(x,cal_func(x,*rguess),'g-')
    plt.errorbar(x,ry,yerr=rye,color='#000000',fmt='.')
    plt.plot(x,cal_func(x,*rpopt),'r-')
fp.format_plot(plt, 1400, 850, 50, 50, tight=False)
plt.show()

# print('Ms*t mean = '+str(np.mean(mstlist))+' +/- '+str(np.std(mstlist)))
print('h mean = '+str(np.mean(hlist))+' +/- '+str(np.std(hlist)))
print('theta mean = '+str(np.mean(thetalist))+' +/- '+str(np.std(thetalist)))
#print('bz0 mean = '+str(np.mean(bz0list))+' +/- '+str(np.std(bz0list)))
