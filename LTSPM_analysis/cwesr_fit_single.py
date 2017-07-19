# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 07:34:44 2016

@author: alec
"""

from scipy.optimize import curve_fit
from scipy import signal
import numpy as np
#import peakdet
import peakutils

dgamp = 1.5e4
dgwidth = 4
pdheight = 1500
maxcenshift = 20

def fit_gaussian(x, *params):
	y = np.zeros_like(x)
	c = params[0]
	freq_ctr = params[1]
	freq_split = params[2]
	amp1 = params[3]
	width1 = params[4]
	amp2 = params[5]
	width2 = params[6]
	y = -abs( amp1 * np.exp(-((x-(freq_ctr + freq_split/2))**2)/(2*(width1**2))) +
				amp2 * np.exp(-((x-(freq_ctr - freq_split/2))**2)/(2*(width2**2))) )
	y=y+c
	return y

def fit_lorentzian(x, *params):
	y = np.zeros_like(x)
	c = params[0]
	for i in range(1, len(params), 3):
		ctr = params[i]
		amp = params[i+1]
		wid = params[i+2]
		#y = y + (amp * (np.exp( -((x - ctr-2.3)/wid)**2)+np.exp( -((x - ctr)/wid)**2)+np.exp( -((x - ctr+2.3)/wid)**2)) + c
		y = y + (-abs(amp * (wid/2)**2)/((x-ctr)**2+(wid/2)**2))
	y=y+c
	return y

def cwesr_fit(x, y, filenum=0, gauss=True, gamp=dgamp, gwidth=dgwidth, gctr=2870, d_gsplit=20):

	dlen = len(y)
	b, a = signal.butter(1, 0.5, btype='lowpass')
	yfilt = signal.filtfilt(b, a, y)
	indexes = peakutils.indexes(np.max(yfilt)-yfilt, thres=0.45, min_dist=2)
	sm1=[]
	sm2=[]
	mintab = np.transpose([x[indexes], y[indexes]])

	for i in range(0,len(mintab)):
		if mintab[i][0] < 2872:
			sm1.append(mintab[i])
		else:
			sm2.append(mintab[i])
	sm1s=sorted(sm1, key=lambda x: x[1])
	sm2s=sorted(sm2, key=lambda x: x[1])

	dmax = np.max(y)
	dmin = np.min(y)
	amp = dmax-dmin

	if len(sm1s) >= 1 and len(sm2s) >= 1:
		fc1 = sm1s[0][0]
		fc2 = sm2s[0][0]
		gsplit = np.abs(fc2-fc1)
		gamp = y[0]-sm1s[0][1]
	elif len(sm1s) == 0 and len(sm2s) >= 1:
		gsplit=0
		gamp = (y[0]-sm2s[0][1])/2
	elif len(sm1s) >= 1 and len(sm2s) == 0:
		gsplit=0
		gamp = (y[0]-sm1s[0][1])/2
	else:
		gsplit = d_gsplit
		gamp = y[0]-dmin

	lbounds2 = [0,gctr,0,amp/3,4,amp/3,4]
	ubounds2 = [3e5,2875,200,2*amp,15,2*amp,15]
	guess = [y[0], gctr, gsplit, gamp, gwidth, gamp, gwidth]

	try:
		if (gauss):
			popt, pcov = curve_fit(fit_gaussian, x, y, p0=guess, bounds=(lbounds2,ubounds2))
		else:
			popt, pcov = curve_fit(fit_lorentzian, x, y, p0=guess, bounds=(lbounds2,ubounds2))
	except:
		popt = [0, 0, 1e3, 1, 1, 1, 1]
		pcov = np.zeros((7,7))
		print('fit fail on file '+str(filenum))

	if (gauss):
		fit = fit_gaussian(x, *popt)
		fitg = fit_gaussian(x, *guess)
	else:
		fit = fit_lorentzian(x, *popt)
		fitg = fit_lorentzian(x, *guess)

	return popt, pcov, fit, fitg, np.transpose(mintab)
