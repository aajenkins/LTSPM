# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np
import scipy.fftpack as fft

#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

def vector_reconstruction(data, theta, phi, height, scansize, kcutoff=1):
	pi = np.pi
	fdata = fft.fft2(data)
	fdata = fft.fftshift(fdata)

	dlen = len(fdata)
	hlen = int(np.floor(dlen/2))

	hxf = np.zeros_like(fdata)
	hyf = np.zeros_like(fdata)
	hzf = np.zeros_like(fdata)
	meffk = np.zeros_like(fdata)
	k = 0
	kmax = 2*pi*dlen/scansize

	for j in range(0,dlen):
	    ky = 2*pi*(j-hlen)/scansize
	    for i in range(0,dlen):
	        kx = 2*pi*(i-hlen)/scansize
	        k = np.sqrt(kx**2 + ky**2)
	        if (i==hlen and j==hlen):
	            hzf[j,i] = fdata[j,i]/np.sin(theta)
	            meffk[j,i] = 0
	        else:
	            hzf[j,i] = fdata[j,i]/(np.cos(theta)*(1-
	            (1j/k)*np.tan(theta)*(kx*np.cos(phi) + ky*np.sin(phi))))
	            hxf[j, i] = -1j*(kx/k)*hzf[j,i]
	            hyf[j, i] = -1j*(ky/k)*hzf[j,i]
	            if (k<kmax*kcutoff):
	                meffk[j,i] = (1/(k))*np.exp(height*k)*hzf[j,i]

	bzdata = np.real(fft.ifft2(fft.ifftshift(hzf)))
	bxdata = np.real(fft.ifft2(fft.ifftshift(hxf)))
	bydata = np.real(fft.ifft2(fft.ifftshift(hyf)))
	meffdata = np.real(fft.ifft2(fft.ifftshift(meffk)))

	return bxdata, bydata, bzdata, meffdata