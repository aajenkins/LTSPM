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

dgamp = 1.4e4
dgwidth = 4

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
	freq_ctr = params[1]
	freq_split = np.abs(params[2])
	amp1 = params[3]
	width1 = params[4]
	amp2 = params[5]
	width2 = params[6]
	y = -abs( (amp1 * (width1)**2)/((x-(freq_ctr + freq_split/2))**2+(width1)**2) +
			(amp2 * (width2)**2)/((x-(freq_ctr - freq_split/2))**2+(width2)**2) )
	y=y+c
	return y

def cwesr_fit(x, y, filenum=0, gauss=False, gamp=dgamp, gwidth=dgwidth, gctr=2870, d_gsplit=20,
				min_width=4, max_width=15, max_counts=3e5, max_ctr=2875, max_splitting=200):

	dlen = len(y)
	b, a = signal.butter(1, 0.4, btype='lowpass')
	yfilt = signal.filtfilt(b, a, y)
	indexes = peakutils.indexes(np.max(yfilt)-yfilt, thres=0.45, min_dist=2)
	sm1=[]
	sm2=[]
	mintab = np.transpose([x[indexes], y[indexes]])

	for i in range(0,len(mintab)):
		if mintab[i][0] < 2873:
			sm1.append(mintab[i])
		else:
			sm2.append(mintab[i])
	sm1s=sorted(sm1, key=lambda x: x[1])
	sm2s=sorted(sm2, key=lambda x: x[1])

	dmax = np.max(y)
	dmin = np.min(y)
	amp = dmax-dmin

	ystart = np.mean(np.append(y[0:3], y[-3:]))

	if len(sm1s) >= 1 and len(sm2s) >= 1:
		fc1 = sm1s[0][0]
		fc2 = sm2s[0][0]
		gsplit = np.abs(fc2-fc1)
		# gctr = (fc2+fc1)/2
		gamp = ystart-sm1s[0][1]
	elif len(sm1s) == 0 and len(sm2s) == 1:
		gsplit=0
		gamp = (ystart-sm2s[0][1])/2
		gctr = sm2s[0][0]
	# elif len(sm1s) == 0 and len(sm2s) >= 2:
	# 	fc1 = sm2s[0][0]
	# 	fc2 = sm2s[1][0]
	# 	gsplit = np.abs(fc2-fc1)
	# 	gamp = (ystart-sm2s[0][1])
	# 	gctr = (fc2+fc1)/2
	elif len(sm1s) == 1 and len(sm2s) == 0:
		gsplit=0
		gamp = (ystart-sm1s[0][1])/2
		gctr = sm1s[0][0]
	# elif len(sm1s) >= 2 and len(sm2s) == 0:
	# 	fc1 = sm1s[0][0]
	# 	fc2 = sm1s[1][0]
	# 	gsplit = np.abs(fc2-fc1)
	# 	gamp = (ystart-sm1s[0][1])
	# 	gctr = (fc2+fc1)/2
	else:
		gsplit = d_gsplit
		gamp = ystart-dmin

	lbounds2 = [0,2865,-max_splitting,amp/8,min_width,amp/8,min_width]
	ubounds2 = [max_counts,max_ctr,max_splitting,4*amp,max_width,4*amp,max_width]
	guess = [ystart, gctr, gsplit, gamp, gwidth, gamp, gwidth]

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

	return popt, pcov, fit, fitg
