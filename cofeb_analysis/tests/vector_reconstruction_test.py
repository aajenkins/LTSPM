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
import scipy.fftpack as fft

import load_scan as lscan
import vector_reconstruction as vr
import fourier_image as fi
import format_plot as fp

pi = np.pi

xres = 250
yres = 250
zfield = 9.5
scanL = 2.5e-4


path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
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

testpath = '/Users/alec/UCSB/mathematica/CoFeB-MgO/skyrmion_shape/tests/'
bzdata = np.loadtxt(testpath+'b_z.txt', delimiter=',')

#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

hzf = fft.fftshift(fft.fft2(bzdata))

dlen = len(hzf)
hlen = int(np.floor(dlen/2))

Vk = np.zeros_like(hzf)
k = 0
kmax = 2*pi*dlen/scanL

for j in range(0,dlen):
	ky = 2*pi*(j-hlen)/scanL
	for i in range(0,dlen):
		kx = 2*pi*(i-hlen)/scanL
		k = np.sqrt(kx**2 + ky**2)
		if (i==hlen and j==hlen):
			Vk[j,i] = 0
		else:
			Vk[j, i] = -hzf[j,i]/(k**2)

Vdata = np.real(fft.ifft2(fft.ifftshift(Vk)))

np.savetxt(testpath+'Vdata.txt', Vdata, delimiter=',')


#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = ax1.imshow(bzdata, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 550, 50)

fig1, ax1 = plt.subplots()
im1 = ax1.imshow(Vdata, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

plt.show()
