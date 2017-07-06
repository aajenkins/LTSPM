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
phinum = 16
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
xy0s = (r0s*np.cos(phis), r0s*np.sin(phis))


np.savetxt(path+'stray_field_sim/'+'radiusphi_'+filespec+'.txt', (angles, r0s), delimiter = ',')

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

# mdataflat = np.array([])
# for i in range(0,len(mdata)):
#     mdataflat = np.append(mdataflat, mdata[i])
#
# fig, ax = plt.subplots()
# plt.hist(mdataflat, bins=60)
# fp.format_plot(plt, 500, 500, 450, 250)


# fig3, ax3 = plt.subplots(nrows=phinum, sharex=True, sharey=True)
# for i in range(0,phinum):
#     ax3[i].plot(mphi[i],'b.')
#     ax3[i].plot(xf, fit_tanh(xf, *guesses[i]), 'g')
#     ax3[i].plot(xf, fits[i], 'r')
#     ax3[i].get_yaxis().set_visible(False)
#
# fig3.subplots_adjust(hspace=0)
# plt.setp([a.get_xticklabels() for a in fig3.axes[:-1]], visible=False)
# fp.format_plot(plt, 400, 800, 850, 50, tight=False)

# fig3, ax3 = plt.subplots(figsize=(4,3.5))
# ax3.plot(xf*0.06, mphi[8],'b.')
# ax3.plot(xf*0.06, fits[8], 'r')
# # plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# # ax3.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# # ax3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax3.set_xticks(np.linspace(0.0,0.9,5))
# # ax1.set_yticks(np.linspace(0.88,1.0,4))
# fp.format_plot(plt, 400, 400, 450, 50)
# pylab.savefig('/Users/alec/UCSB/scan_images/phicut_'+str(scannum)+filespec+'.png', format='png')


# fig4, ax4 = plt.subplots()
# im4 = plt.imshow(mdataint, cmap='jet', interpolation='nearest')
# plt.colorbar(im4, fraction=0.046, pad=0.04)
# fp.format_plot(plt, 400, 400, 50, 450)

# fig4, ax4 = plt.subplots()
# im4 = plt.imshow(bzdataError, cmap='jet', interpolation='nearest')
# plt.colorbar(im4, fraction=0.046, pad=0.04)
# fp.format_plot(plt, 400, 400, 50, 450)


# fig5, ax5 = plt.subplots()
# im5 = plt.imshow(bzdata, cmap='gray', interpolation='nearest')
# plt.colorbar(im5, fraction=0.046, pad=0.04)
# fp.format_plot(plt, 400, 400, 450, 450)


# fig6, ax6 = plt.subplots()
# x = np.linspace(0,2*pi,phinum)
# plt.plot(x, widths)
# fp.format_plot(plt, 400, 400, 450, 450)

plt.show()
