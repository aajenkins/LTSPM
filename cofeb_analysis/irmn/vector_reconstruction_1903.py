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
scanL = 0.5*5000

savepath = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
data = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt', xres, yres, 15)
#misc.imsave('/Users/alec/UCSB/scan_images/full-field/ff1760.png', data[0])
ffmask = ndimage.imread('/Users/alec/UCSB/scan_images/full-field/ff1903mask.png', flatten=True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

path = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
filespec = 'Msnotfixed'
cal_params = np.loadtxt(path+'cal_parameters_'+filespec+'.txt', delimiter=',')

Ms = cal_params[0]
t = cal_params[1]
theta = cal_params[2]
phi = cal_params[3]
height = cal_params[4]
heighterr = cal_params[5]

#---------------- FIT FUNCTIONS ----------------------------------
#-----------------------------------------------------------------

def fit_arctan(x, *params):
    y = np.zeros_like(x)
    c = params[0]
    a = params[1]
    x0 = params[2]
    wid = params[3]

    y = c+(a/pi)*np.arctan((x-x0)/wid)
    return y

#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

datas = np.multiply(ffmask,data[0])
datas0 = np.add(datas,-np.cos(theta)*zfield)
datas0filter = signal.wiener(datas0, 5)

recon_data = vr.vector_reconstruction(datas0filter, theta, phi, height, scanL, kcutoff=0.3)

bxdata = recon_data[0]
bydata = recon_data[1]
bzdata = recon_data[2]
mdata = recon_data[3]

mdataint = ndimage.interpolation.zoom(mdata, 2, order=1)

np.savetxt(savepath+'bz_'+str(scannum)+filespec+'.txt', bzdata, delimiter=',')

#---------------- LINECUTS ---------------------------------------
#-----------------------------------------------------------------

x0, y0 = 21.5, 28
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
	popt, pcov = curve_fit(fit_arctan, xf, mphi[i], p0=guesses[i])
	fits[i] = fit_arctan(xf, *popt)
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
    bpopt, bpcov = curve_fit(fit_arctan, bxf, bzphi[i], p0=bguesses[i])
    bfits[i] = fit_arctan(xf, *bpopt)
    bwidths[i] = np.abs(bpopt[3])

print(np.mean(widths))
print(np.std(widths))

#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(datas0, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(datas0filter, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mdata, cmap='jet', interpolation='nearest')
plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, use_gridspec=True)

for i in range(0,phinum):
    phi = i*2*pi/phinum
    x1, y1 = x0-lclen*np.cos(phi), y0-lclen*np.sin(phi)
    plt.plot([x0, x1], [y0, y1], 'k-')
plt.axis('image')
fp.format_plot(plt, 400, 400, 450, 50)


fig, axes = plt.subplots(nrows=phinum, sharex=True, sharey=True)
for i in range(0,phinum):
	axes[i].plot(mphi[i],'b.')
	axes[i].plot(xf, fit_arctan(xf, *guesses[i]), 'g')
	axes[i].plot(xf, fits[i], 'r')
	axes[i].get_yaxis().set_visible(False)

fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
#fig.tight_layout()

fp.format_plot(plt, 400, 900, 900, 50, tight=False)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mdataint, cmap='jet', interpolation='nearest')
plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 450)


fig1, ax1 = plt.subplots()
plt.imshow(bzdata, cmap='gray', interpolation='nearest')
plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 450, 450)
#for i in range(0,phinum):
#    phi = i*2*pi/phinum
#    x1, y1 = x0-lclen*np.cos(phi), y0-lclen*np.sin(phi)
#    plt.plot([x0, x1], [y0, y1], 'b-')
#plt.axis('image')

# fp.format_plot(plt, 50, 450)
#pylab.savefig('/Users/alec/UCSB/scan_images/tacofeb_reconstruction/recon_bz_'+str(scannum)+'.pdf')

#fig, axes = plt.subplots(nrows=4, sharex=True, sharey=True, figsize=(5,9))
#for i in range(0,4):
#    axes[i].plot(bzphi[i],'b.')
##    axes[i].plot(xf,fits[i], 'r')
#    axes[i].get_yaxis().set_visible(False)
#
#fig.subplots_adjust(hspace=0)
#plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
##fig.tight_layout()
#
#fig = plt.gcf().canvas.manager.window
#geom = fig.geometry()
#x,y,dx,dy = geom.getRect()
#fig.setGeometry(900,50,dx, dy)
#fig.raise_()
#
fig1, ax1 = plt.subplots()
x = np.linspace(0,2*pi,phinum)
plt.plot(x, r0s)
fp.format_plot(plt, 400, 400, 850, 450)

plt.show()
