# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np

#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

def vector_reconstruction(data, dataError, theta, thetaError, phi, height, scansize, t, kcutoff=1):
	pi = np.pi
	fdata = np.fft.fft2(data)
	fdata = np.fft.fftshift(fdata)

	dlen = len(fdata)
	hlen = int(np.floor(dlen/2))

	hxf = np.zeros_like(fdata)
	hyf = np.zeros_like(fdata)
	hzf = np.zeros_like(fdata)
	mzBlochk = np.zeros_like(fdata)
	Vk = np.zeros_like(fdata)
	k = 0
	kmax = 2*pi*dlen/scansize

	for j in range(0,dlen):
		ky = 2*pi*(j-hlen)/scansize
		for i in range(0,dlen):
			kx = 2*pi*(i-hlen)/scansize
			k = np.sqrt(kx**2 + ky**2)
			if (i==hlen and j==hlen):
				hzf[j,i] = 0#fdata[j,i]/np.sin(theta)
				mzBlochk[j,i] = 0
				Vk[j,i] = 0
			else:
				if (k<kmax*kcutoff):
					hzf[j,i] = fdata[j,i]/( np.cos(theta) * (1-
						(1j/k)*np.tan(theta)*(kx*np.cos(phi) + ky*np.sin(phi))) )
					hxf[j, i] = -1j*(kx/k)*hzf[j,i]
					hyf[j, i] = -1j*(ky/k)*hzf[j,i]
					Vk[j, i] = -hzf[j,i]/(k**2)
					mzBlochk[j,i] = 2*np.exp(height*k)*hzf[j,i]/(1-np.exp(-k*t))

	bzdata = np.real(np.fft.ifft2(np.fft.ifftshift(hzf)))
	bxdata = np.real(np.fft.ifft2(np.fft.ifftshift(hxf)))
	bydata = np.real(np.fft.ifft2(np.fft.ifftshift(hyf)))

	bzdataError = np.zeros_like(bzdata)
	for j in range(0,dlen):
		for i in range(0,dlen):
			bzthetaError = ( (np.sin(theta)/(np.cos(theta)**2)*(data[j,i] - bxdata[j,i]*np.sin(theta)*np.cos(phi) - bydata[j,i]*np.sin(theta)*np.sin(phi)))
							+ (1/np.cos(theta))*(-bxdata[j,i]*np.cos(theta)*np.cos(phi) - bydata[j,i]*np.cos(theta)*np.sin(phi)) )*thetaError
			bzbnvError = np.cos(theta)*dataError[j,i]
			bzdataError[j,i] = np.sqrt((bzthetaError)**2 + (bzbnvError)**2)

	mzBlochdata = np.real(np.fft.ifft2(np.fft.ifftshift(mzBlochk)))
	Vdata = np.real(np.fft.ifft2(np.fft.ifftshift(Vk)))

	return bxdata, bydata, bzdata, mzBlochdata, Vdata, bzdataError, mzBlochk

def vector_reconstruction_1D(data, theta, phi):
	# reconstruct two field components for field uniform along one xaxis_scale

	pi = np.pi
	bnvf = np.fft.fft(data)
	bnvf = np.fft.fftshift(bnvf)

	dlen = len(bnvf)
	hlen = int(np.floor(dlen/2))

	bxf = np.zeros_like(bnvf)
	bzf = np.zeros_like(bnvf)

	for i in range(0,dlen):
		kx = 2*pi*(i-hlen)/dlen
		if (i==hlen):
			bxf[i] = 0
			bzf[i] = 0
		else:
			bxf[i] = bnvf[i]/( np.sin(theta) * np.cos(phi) - 1j * np.sign(kx) * np.cos(theta) )
			bzf[i] = -1j * np.sign(kx) * bxf[i]

	bxdata = np.real(np.fft.ifft(np.fft.ifftshift(bxf)))
	bzdata = np.real(np.fft.ifft(np.fft.ifftshift(bzf)))

	return bxdata, bzdata
