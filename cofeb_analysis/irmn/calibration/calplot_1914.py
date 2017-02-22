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

matplotlib.rc('font', **font)

def cal_func(x, *args):
    bnv = np.zeros_like(x)
    h = args[0]
    x0 = args[1]
    theta = args[2]*pi/180

    bx = (2e-3)*mstnm*(h/(h**2+(x-x0)**2))
    bz = -(2e-3)*mstnm*((x-x0)/(h**2+(x-x0)**2))
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

ploth = 1
fig1, ax1 = plt.subplots()
gs = gridspec.GridSpec(5, 8)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)

guess = [120, 1500, 55]
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
	hlist[2*j] = popt[0]
	hlist[2*j+1] = rpopt[0]
	thetalist[2*j] = popt[2]
	thetalist[2*j+1] = rpopt[2]
	#    bz0list[2*j] = popt[4]
	#    bz0list[2*j+1] = rpopt[4]
	#

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
