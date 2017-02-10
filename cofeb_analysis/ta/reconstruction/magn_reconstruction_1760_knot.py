# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np
from scipy import ndimage
from scipy import signal
import load_scan as lscan
import vector_reconstruction as vr
import fourier_image as fi
import m_reconstruction_Dovzhenko as mr

pi = np.pi

scannum = 1760
xres = 50
yres = 50
zfield = 9.5
scanL = 0.6*5000

data = lscan.load_ff('fitdata_1760.txt',xres,yres,15)
ffmask = ndimage.imread('ff_1760mask.png',flatten=True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

filespec = 'Msnotfixed'
cal_params = np.loadtxt('cal_parameters_'+filespec+'.txt', delimiter=',')

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
mzdataintnorm = signal.wiener(mzdataintnorm)

phi = mr.m_reconstruction_Dovzhenko(mzdataintnorm, 400, 2*pi/180)

for i in range(0,len(phi)):
    for j in range(0,len(phi)):
        phi[j][i] = (np.add(np.multiply(phi[j][i],180/pi),180))%36

np.savetxt('phi.txt', phi, delimiter=',')
