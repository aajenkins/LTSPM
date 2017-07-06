# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import misc
from scipy import signal
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pylab as pylab
from matplotlib.ticker import FormatStrFormatter
import json
import shapely.geometry as geom

import load_scan as lscan
import vector_reconstruction as vr
import fourier_image as fi
import format_plot as fp

pi = np.pi

scannum = 1760
xres = 50
yres = 50
zfield = 9.5
scanL = 0.6*(5e-6)


path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'+str(scannum)+'/'
filespec = 'Msfixed'
cal_params_path = path+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

Ms = cal_params['Ms']
t = cal_params['t']
phi = cal_params['phi']
theta = cal_params['theta']
thetaError = cal_params['thetaError']
height = cal_params['height']
heightError = cal_params['heightError']


data = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,20)
misc.imsave('/Users/alec/UCSB/scan_images/full-field/ff1760.png', data[0])
ffmask = ndimage.imread('/Users/alec/UCSB/scan_images/full-field/ff1760mask.png',flatten = True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

#---------------- FIT FUNCTIONS ----------------------------------
#-----------------------------------------------------------------

def fit_tanh(x, *params):
    y = np.zeros_like(x)
    c = params[0]
    a = params[1]
    x0 = params[2]
    wid = params[3]

    y = c+(a/2)*np.tanh((x-x0)/wid)
    return y


def get_wall_distance(x, y, phis_int_list, r0s_int):
    x0 = r0s_int(phis_int_list)*np.cos(phis_int_list)
    y0 = r0s_int(phis_int_list)*np.sin(phis_int_list)
    coarsedist = np.sqrt((x0-x)**2 + (y0-y)**2)
    minargs = np.argsort(coarsedist)[:2]
    phi0 = phis_int_list[minargs[0]]
    phi1 = phis_int_list[minargs[1]]
    dist0 = coarsedist[minargs[0]]
    dist1 = coarsedist[minargs[1]]
    minargs = phis_int_list[minargs]

    for i in range(0,4):
        midphi = (phi0+phi1)/2
        middist = np.sqrt((r0s_int(midphi)*np.cos(midphi)-x)**2 + (r0s_int(midphi)*np.sin(midphi)-y)**2)
        if (middist < dist1):
            dist1 = middist
            phi1 = midphi
        if (dist1 < dist0):
            tempdist = dist0
            tempphi = phi0
            dist0 = dist1
            phi0 = phi1
            dist1 = tempdist
            phi1 = tempphi

    if (r0s_int(phi0) > np.sqrt(x**2 + y**2)):
        dist0 = -dist0
    return dist0


#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

datas = np.multiply(ffmask,data[0])
datas0 = np.add(datas,-np.cos(theta)*zfield)
datas0filter = signal.wiener(datas0, 3)
datas0filter = signal.medfilt(datas0filter, 5)

wimdata = fi.window_image(datas0)
wimdatafilter = fi.window_image(datas0filter, power=1/2)

recon_data = vr.vector_reconstruction(wimdata, data[1], theta, thetaError, phi, height, scanL, kcutoff=1)

bxdata = recon_data[0]
bydata = recon_data[1]
bzdata = recon_data[2]+zfield
mdata = recon_data[3]
Vdata = recon_data[4]
bzdataError = recon_data[5]

mdataint = ndimage.interpolation.zoom(mdata, 2, order=1)

np.savetxt(path+'bz_'+str(scannum)+'_'+filespec+'.txt', bzdata, delimiter=',')
np.savetxt(path+'bzError_'+str(scannum)+'_'+filespec+'.txt', bzdataError, delimiter=',')
np.savetxt(path+'bx_'+str(scannum)+'_'+filespec+'.txt', bxdata, delimiter=',')
np.savetxt(path+'by_'+str(scannum)+'_'+filespec+'.txt', bydata, delimiter=',')
np.savetxt(path+'mdata_'+str(scannum)+'_'+filespec+'.txt', mdata, delimiter=',')
np.savetxt(path+'V_'+str(scannum)+'_'+filespec+'.txt', Vdata, delimiter=',')

#---------------- LINECUTS ---------------------------------------
#-----------------------------------------------------------------

x0, y0 = 24, 29
phinum = 30
lcnum = 15
lclen = 15
mphi = np.zeros((phinum,lcnum))
for i in range(0,phinum):
    phi = i*2*pi/phinum
    x1, y1 = x0+lclen*np.cos(phi), y0+lclen*np.sin(phi)
    x, y = np.linspace(x0, x1, lcnum), np.linspace(y0, y1, lcnum)
    mphi[i] = ndimage.map_coordinates(np.transpose(mdata), np.vstack((x,y)), order=1)

#---------------- LINE FITS --------------------------------------
#-----------------------------------------------------------------

xf = np.arange(0, lcnum)
fits = np.zeros((phinum, lcnum))
guesses = np.zeros((phinum, 4))
widths = np.zeros(phinum)
r0s = np.zeros(phinum)
angles = np.linspace(0, 2*pi, phinum+1)

for i in range (0,phinum):
    y = mphi[i]
    guesses[i] = [(y[-1]+y[0])/2,y[-1]-y[0],6,1]
    popt, pcov = curve_fit(fit_tanh, xf, mphi[i], p0=guesses[i])
    fits[i] = fit_tanh(xf, *popt)
    widths[i] = np.abs(popt[3])
    r0s[i] = popt[2]*scanL/xres

r0s = np.append(r0s,r0s[0])
phis = np.linspace(0, 2*np.pi, phinum+1)
x0, y0 = (r0s*np.cos(phis), r0s*np.sin(phis))
tck, u = interpolate.splprep([x0, y0], s=0, k=3, per=True)
t = np.zeros(x0.shape)
t[1:] = np.sqrt((x0[1:] - x0[:-1])**2 + (y0[1:] - y0[:-1])**2)
t = np.cumsum(t)
t /= t[-1]
tnew = np.linspace(0, 1, 50)
xnew,ynew = interpolate.splev(tnew, tck)
dw_grads = np.zeros((2, len(xnew)-1))
theta_grad = np.zeros((len(xnew)-1))
theta = np.arctan2(ynew[:-1], xnew[:-1])

dw_grads[0] = -(ynew[1:]-ynew[0:-1])
dw_grads[1] = (xnew[1:]-xnew[0:-1])
gradlens = np.sqrt(dw_grads[0]**2 + dw_grads[1]**2)
dw_grads[0] = np.divide(dw_grads[0], gradlens)
dw_grads[1] = np.divide(dw_grads[1], gradlens)
theta_grad = np.arctan(dw_grads[1], dw_grads[0])

theta_grad_int = interpolate.interp1d(theta, theta_grad, bounds_error=False,
                                      fill_value=0., kind='cubic')
theta_int_list = np.linspace(min(theta),max(theta),100)
r0s_int = interpolate.interp1d(phis, r0s, bounds_error=False,
                                      fill_value=0., kind='cubic')
phis_int_list = np.linspace(min(phis),max(phis),100)
x0_int, y0_int = (r0s_int(phis_int_list)*np.cos(phis_int_list), r0s_int(phis_int_list)*np.sin(phis_int_list))

gridnum = 200
xm_grid, ym_grid = np.linspace(-scanL/2, scanL/2, gridnum), np.linspace(-scanL/2, scanL/2, gridnum)

domain_wall = geom.LineString(np.transpose([x0_int, y0_int]))

grid_dist = np.zeros((gridnum, gridnum))
for j in range(gridnum):
    for i in range(gridnum):
        point = geom.Point(xm_grid[i], ym_grid[j])
        grid_dist[j, i] = domain_wall.distance(point)

# BLOCH
dw = 30.0e-9
mz = np.tanh(grid_dist/dw)


#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

fig2, ax2 = plt.subplots(figsize=(5,5))
im2 = plt.imshow(mdata, cmap='gray', interpolation='nearest')
plt.colorbar(im2, fraction=0.046, pad=0.04, use_gridspec=True)
ax2.set_axis_off()
fp.format_plot(plt, 350, 350, 450, 50)

fig3, ax3 = plt.subplots()
plt.plot(theta, theta_grad, 'g.')
plt.plot(theta_int_list, theta_grad_int(theta_int_list), 'r-')
fp.format_plot(plt, 350, 350, 50, 50)

fig3, ax3 = plt.subplots()
plt.plot(phis, r0s, 'g.')
plt.plot(phis_int_list, r0s_int(phis_int_list), 'r-')
fp.format_plot(plt, 350, 350, 50, 450)

fig4, ax4 = plt.subplots()
plt.plot(x0, y0, 'g.')
plt.plot(x0_int, y0_int, 'r-')
fp.format_plot(plt, 350, 350, 450, 450)

fig4, ax4 = plt.subplots()
im = plt.imshow(mz)
plt.colorbar(im, fraction=0.046, pad=0.04, use_gridspec=True)
fp.format_plot(plt, 350, 350, 450, 450)


plt.show()
