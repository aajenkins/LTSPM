# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np
from scipy import ndimage
from scipy import misc
from scipy import signal
import load_scan as lscan
import vector_reconstruction as vr
import fourier_image as fi
import format_plot as fp
import m_reconstruction_Dovzhenko as mr

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

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
filespec = 'Msnotfixed'
cal_params = np.loadtxt(path+'cal_parameters_'+filespec+'.txt', delimiter=',')

theta = cal_params[2]
phi = cal_params[3]
height = cal_params[4]

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

mzdatafilt= signal.medfilt(meffdata, 7)
mzdataint = ndimage.interpolation.zoom(mzdatafilt, 2, order=1)

minmz = np.min(mzdataint)
maxmz = np.max(mzdataint)
mzdataintnorm = np.multiply(np.add(mzdataint,-(maxmz+minmz)/2),1.999999/(maxmz-minmz))

np.savetxt(path+'mzdata_norm.txt',mzdataintnorm,delimiter=',')

phi = mr.m_reconstruction_Dovzhenko(mzdataintnorm, 50, 0.005, 3*pi/180)

np.savetxt(path+'phi.txt',phi,delimiter=',')
