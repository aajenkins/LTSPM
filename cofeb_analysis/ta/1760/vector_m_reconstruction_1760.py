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
numxys = 10
scales = np.linspace(3.0, 0.2, numxys)
meanrs = np.zeros(numxys)
res = np.zeros((numxys, len(xnew)))
coffset = 0.4
xys = np.zeros((numxys, 2, len(xnew)))
dw_grads = np.zeros((numxys, 2, len(xnew)-1))
theta = np.zeros((numxys, len(xnew)-1))


for i in range(numxys):
    xys[i, 0] = xnew*scales[i]
    xys[i, 1] = ynew*scales[i]
    dw_grads[i, 0] = -(xys[i, 1,1:]-xys[i, 1,0:-1])
    dw_grads[i, 1] = (xys[i, 0,1:]-xys[i, 0,0:-1])
    gradlens = np.sqrt(dw_grads[i, 0]**2 + dw_grads[i, 1]**2)
    dw_grads[i, 0] = np.divide(dw_grads[i, 0], gradlens)
    dw_grads[i, 1] = np.divide(dw_grads[i, 1], gradlens)
    # theta[i] = np.arctan2(dw_grads[i, 1], dw_grads[i, 0])
    # thetas[i, 2] = np.arctan2(dw_grads[i, 1], dw_grads[i, 0])
    # meanrs[i] = np.mean(np.sqrt(xys[i, 0]**2 + xys[i, 1]**2))
    # res[i] = meanrs[i] - np.sqrt(xys[i, 0]**2 + xys[i, 1]**2)
    # thetas[i] = np.arctan2(xys[i, 1], xys[i, 0])
    # xys[i, 0] = meanrs[i]*np.cos(thetas[i])*(1-scales[i]) + xys[i, 0]*(scales[i])
    # xys[i, 1] = meanrs[i]*np.sin(thetas[i])*(1-scales[i]) + xys[i, 1]*(scales[i])

# thetas = [xys[:, 0, 0:-1].flatten(), xys[:, 1, 0:-1].flatten(),
            #    theta.flatten()]

points = np.transpose([xys[:, 0, 0:-1].flatten(), xys[:, 1, 0:-1].flatten()])
dw_gradx = dw_grads[:, 0].flatten()
dw_grady = dw_grads[:, 1].flatten()

#np.savetxt(path+'stray_field_sim/'+'radiusphi_'+filespec+'.txt', (angles, r0s), delimiter = ',')

gridx, gridy = np.mgrid[-4.5e-7:4.5e-7:200j, -4.5e-7:4.5e-7:200j]
dw_gradx_int = interpolate.griddata(points, dw_gradx, (gridx, gridy), method='linear')
dw_grady_int = interpolate.griddata(points, dw_grady, (gridx, gridy), method='linear')
thetaint = np.arctan2(dw_grady_int, dw_gradx_int)
thetaint = np.mod(thetaint, 2*pi)

#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = ax1.imshow(wimdatafilter, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

fig2, ax2 = plt.subplots(figsize=(5,5))
im2 = plt.imshow(mdata, cmap='jet', interpolation='nearest')
plt.colorbar(im2, fraction=0.046, pad=0.04, use_gridspec=True)
ax2.set_axis_off()
# for i in range(0,phinum):
#     phi = i*2*pi/phinum
#     x1, y1 = x0+lclen*np.cos(phi), y0+lclen*np.sin(phi)
#     plt.plot([x0, x1], [y0, y1], 'k-')
# plt.axis('image')
fp.format_plot(plt, 350, 350, 450, 50)
# pylab.savefig('/Users/alec/UCSB/scan_images/mdata_'+str(scannum)+filespec+'.png', format='png')

fig3, ax3 = plt.subplots(figsize=(5,5))
#plt.plot(y0, x0)
#plt.plot(ynew, xnew, color='r')
for i in range(numxys):
    plt.plot(xys[i, 0], xys[i, 1], color='r')
ax3.set_aspect('equal', 'datalim')
fp.format_plot(plt, 350, 350, 50, 450)

fig5, ax5 = plt.subplots(figsize=(5,5))
plt.imshow(dw_gradx_int)
fp.format_plot(plt, 350, 350, 450, 450)

fig5, ax5 = plt.subplots(figsize=(5,5))
plt.imshow(dw_grady_int)
fp.format_plot(plt, 350, 350, 50, 450)


plt.show()
