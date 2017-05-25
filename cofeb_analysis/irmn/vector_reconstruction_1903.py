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
from scipy.optimize import curve_fit
import scipy.fftpack as fft
import json
import matplotlib.pylab as pylab
import vector_reconstruction as vr
import load_scan as lscan
import fourier_image as fi
import format_plot as fp

pi = np.pi

scannum = 1903
#scanbacknum = 1740
xres = 50
yres = 50
zfield = 0
scanL = 0.5*5e-6

savepath = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
data = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt', xres, yres, 15)
#misc.imsave('/Users/alec/UCSB/scan_images/full-field/ff1760.png', data[0])
ffmask = ndimage.imread('/Users/alec/UCSB/scan_images/full-field/ff1903mask.png', flatten=True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

path = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
filespec = 'Msfixed'
cal_params_path = path+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

Ms = cal_params['Ms']
t = cal_params['t']
phi = cal_params['phi']
theta = cal_params['theta']
thetaError = cal_params['thetaError']
# height = cal_params['height']
height = 135.0*1e-9
heightError = cal_params['heightError']

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
dlen = len(datas)
coswindow = np.zeros((dlen,dlen))
for j in range(0,dlen):
    for i in range(0,dlen):
        coswindow[j,i] = np.sin(pi*j/(dlen-1)) * np.sin(pi*i/(dlen-1))
# datas0filter = signal.medfilt(datas0, 7)

recon_data = vr.vector_reconstruction(datas0filter, data[1], theta, thetaError, phi, height, scanL, kcutoff=1)

bxdata = recon_data[0]
bydata = recon_data[1]
bzdata = recon_data[2]
mdatak = recon_data[6]*(coswindow**2)
Vdata = recon_data[4]
bzdataError = recon_data[5]

mdata = np.real(fft.ifft2(fft.ifftshift(mdatak)))

mdataint = ndimage.interpolation.zoom(mdata, 2, order=1)

np.savetxt(savepath+'bx_'+str(scannum)+'_'+filespec+'.txt', bxdata, delimiter=',')
np.savetxt(savepath+'by_'+str(scannum)+'_'+filespec+'.txt', bydata, delimiter=',')
np.savetxt(savepath+'bz_'+str(scannum)+'_'+filespec+'.txt', bzdata, delimiter=',')
np.savetxt(path+'bzError_'+str(scannum)+'_'+filespec+'.txt', bzdataError, delimiter=',')
np.savetxt(savepath+'bnv_'+str(scannum)+'_'+filespec+'.txt', data[0], delimiter=',')
np.savetxt(path+'mdata_'+str(scannum)+'_'+filespec+'.txt', mdata, delimiter=',')
np.savetxt(path+'V_'+str(scannum)+'_'+filespec+'.txt', Vdata, delimiter=',')

#---------------- LINECUTS ---------------------------------------
#-----------------------------------------------------------------

x0, y0 = 22, 28
phinum = 16
lcnum = 15
lclen = 15
mphi = np.zeros((phinum,lcnum))
for i in range(0,phinum):
    phi = i*2*pi/phinum
    x1, y1 = x0-lclen*np.cos(phi), y0-lclen*np.sin(phi)
    x, y = np.linspace(x0, x1, lcnum), np.linspace(y0, y1, lcnum)
    mphi[i] = ndimage.map_coordinates(np.transpose(mdata), np.vstack((x,y)), order=1)

x0, y0 = 21.5, 28
phinum = 16
bzlcnum = 15
bzlclen = 15
bzphi = np.zeros((phinum,bzlcnum))
for i in range(0,phinum):
   phi = i*2*pi/phinum
   x1, y1 = x0-lclen*np.cos(phi), y0-lclen*np.sin(phi)
   x, y = np.linspace(x0, x1, bzlcnum), np.linspace(y0, y1, bzlcnum)
   bzphi[i] = ndimage.map_coordinates(np.transpose(bzdata), np.vstack((x,y)), order=1)

#---------------- LINE FITS --------------------------------------
#-----------------------------------------------------------------

xf = np.arange(0,lcnum)
fits = np.zeros((phinum,lcnum))
guesses = np.zeros((phinum, 4))
widths = np.zeros(phinum)
r0s = np.zeros(phinum)
angles = np.linspace(0,2*pi,phinum)

for i in range (0,phinum):
	y = mphi[i]
	guesses[i] = [(y[-1]+y[0])/2,y[-1]-y[0],6,1]
	popt, pcov = curve_fit(fit_tanh, xf, mphi[i], p0=guesses[i])
	fits[i] = fit_tanh(xf, *popt)
	widths[i] = np.abs(popt[3])
	r0s[i] = popt[2]*scanL/xres

np.savetxt('/Users/alec/UCSB/cofeb_analysis_data/irmn/stray_field_sim/radiusphi_'+filespec+'.txt',(angles,r0s), delimiter=',')

bguesses = np.zeros((phinum, 4))
bfits = np.zeros((phinum,bzlcnum))
bwidths = np.zeros(phinum)
bxf = np.arange(0,bzlclen)

for i in range (0,phinum):
    y = bzphi[i]
    bguesses[i] = [(y[-1]+y[0])/2,y[-1]-y[0],6,1]
    bpopt, bpcov = curve_fit(fit_tanh, bxf, bzphi[i], p0=bguesses[i])
    bfits[i] = fit_tanh(xf, *bpopt)
    bwidths[i] = np.abs(bpopt[3])


#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(datas0, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 250, 250, 50, 50)
pylab.savefig('/Users/alec/UCSB/scan_images/datas_'+str(scannum)+filespec+'.png', format='png')

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(bxdata, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 250, 250, 50, 50)
pylab.savefig('/Users/alec/UCSB/scan_images/bx_'+str(scannum)+filespec+'.png', format='png')

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(bydata, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 250, 250, 50, 50)
pylab.savefig('/Users/alec/UCSB/scan_images/by_'+str(scannum)+filespec+'.png', format='png')

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(bzdata, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 250, 250, 50, 50)
pylab.savefig('/Users/alec/UCSB/scan_images/bz_'+str(scannum)+filespec+'.png', format='png')


fig1, ax1 = plt.subplots()
im1 = plt.imshow(datas0filter, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(mdata, cmap='jet', interpolation='nearest')
plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, use_gridspec=True)

for i in range(0,phinum):
    phi = i*2*pi/phinum
    x1, y1 = x0-lclen*np.cos(phi), y0-lclen*np.sin(phi)
    plt.plot([x0, x1], [y0, y1], 'k-')
plt.axis('image')
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 400, 400, 450, 50)
pylab.savefig('/Users/alec/UCSB/scan_images/meff_'+str(scannum)+filespec+'.png', format='png')



fig, axes = plt.subplots(nrows=phinum, sharex=True, sharey=True)
fig.set_size_inches(3, 5)

for i in range(0,phinum):
	axes[i].plot(mphi[i],'b.')
	axes[i].plot(xf, fit_tanh(xf, *guesses[i]), 'g')
	axes[i].plot(xf, fits[i], 'r')
	axes[i].get_yaxis().set_visible(False)

fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
#fig.tight_layout()
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 400, 900, 900, 50, tight=False)
pylab.savefig('/Users/alec/UCSB/scan_images/meff_linecuts_'+str(scannum)+filespec+'.png', format='png')


fig1, ax1 = plt.subplots()
im1 = plt.imshow(mdataint, cmap='jet', interpolation='nearest')
plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 450)


fig1, ax1 = plt.subplots()
plt.imshow(bzdata, cmap='gray', interpolation='nearest')
plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 450, 450)

fig1, ax1 = plt.subplots()
x = np.linspace(0,2*pi,phinum)
plt.plot(x, r0s)
fp.format_plot(plt, 400, 400, 850, 450)

plt.show()
