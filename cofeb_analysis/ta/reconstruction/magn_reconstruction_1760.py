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
from scipy.optimize import curve_fit
import matplotlib.pylab as pylab
import load_scan as lscan
import vector_reconstruction as vr
import fourier_image as fi
import format_plot as fp
import m_reconstruction as mr

pi = np.pi

scannum = 1760
xres = 50
yres = 50
zfield = 9.5
scanL = 0.6*5000

data = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,15)
misc.imsave('/Users/alec/UCSB/scan_images/full-field/ff1760.png', data[0])
ffmask = ndimage.imread('/Users/alec/UCSB/scan_images/full-field/ff_1760mask.png',flatten=True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

path = '/Users/alec/UCSB/LTSPM/cofeb_analysis/ta/'+str(scannum)+'/'
filespec = 'Msnotfixed'
cal_params = np.loadtxt(path+'cal_parameters_'+filespec+'.txt', delimiter=',')

theta = cal_params[2]
phi = cal_params[3]
height = cal_params[4]

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

theta = 55.469*np.pi/180

datas = np.multiply(ffmask,data[0])
datas0 = np.add(datas,-np.cos(theta)*zfield)

wimdata = fi.window_image(datas0)

recon_data = vr.vector_reconstruction(wimdata, theta, phi, height, scanL)

bxdata = recon_data[0]
bydata = recon_data[1]
bzdata = recon_data[2]
meffdata = recon_data[3]

mzdataint = ndimage.interpolation.zoom(meffdata, 2, order=1)
minmz = np.min(mzdataint)
maxmz = np.max(mzdataint)

mzdataintnorm = np.multiply(np.add(mzdataint,-(maxmz+minmz)/2),1.999999/(maxmz-minmz))

phi = mr.m_recon(mzdataintnorm, 20*pi/180, 100)

for i in range(0,len(phi)):
    for j in range(0,len(phi)):
        phi[j][i] = (np.add(np.multiply(phi[j][i],180/pi),180))%360

#---------------- LINECUTS ---------------------------------------
#-----------------------------------------------------------------



#---------------- LINE FITS --------------------------------------
#-----------------------------------------------------------------



#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(datas, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

fig2, ax2 = plt.subplots()
im2 = plt.imshow(mzdataintnorm, cmap='jet', interpolation='nearest')
plt.colorbar(im2, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 450, 50)

fig3, ax3 = plt.subplots()
im3 = plt.imshow(phi, cmap='jet', interpolation='nearest')
plt.colorbar(im3, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 450)

plt.show()
