# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 07:34:44 2016

@author: alec
"""

from scipy.optimize import curve_fit
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
	for i in range(1, len(params), 3):
		ctr = params[i]
		amp = params[i+1]
		wid = params[i+2]
		#y = y + (amp * (np.exp( -((x - ctr-2.3)/wid)**2)+np.exp( -((x - ctr)/wid)**2)+np.exp( -((x - ctr+2.3)/wid)**2)) + c
		y = y + (-abs(amp * np.exp(-((x-ctr)**2)/(2*(wid**2)))))
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

def cwesr_fit(x,y,filenum=0,gauss=True, gamp=dgamp, gwidth=dgwidth, default1=2860, default2=2880):
#    [maxtab, mintab]=peakdet.peakdet(y, y[0]*0.03, x)
	indexes = peakutils.indexes(np.max(y)-y, thres=0.45, min_dist=2)
	sm1=[]
	sm2=[]
	mintab = np.transpose([x[indexes], y[indexes]])

#    if len(indexes) > 0:
#        mintab = np.transpose([x[indexes], y[indexes]])
	for i in range(0,len(mintab)):
		if mintab[i][0] < 2872:
			sm1.append(mintab[i])
		else:
			sm2.append(mintab[i])
    #    mintabs=sorted(mintab, key=lambda x: x[1])
	sm1s=sorted(sm1, key=lambda x: x[1])
	sm2s=sorted(sm2, key=lambda x: x[1])


    #    if len(mintabs) >= 2:
    #        fc1=mintabs[0][0]
    #        fc2 = mintabs[1][0]
    #    elif len(mintabs)==1:
    #        fc1=mintabs[0][0]
    #        fc2=mintabs[0][0]
    #    else:
    #        fc1 = defaultf1
    #        fc2 = defaultf2

	if len(sm1s) >= 1 and len(sm2s) >= 1:
		fc1 = sm1s[0][0]
		fc2 = sm2s[0][0]
		gamp = y[0]-sm1s[0][1]
	elif len(sm1s) == 0 and len(sm2s) >= 1:
		fc1 = sm2s[0][0]
		fc2 = sm2s[0][0]
		gamp = y[0]-sm2s[0][1]
	elif len(sm1s) >= 1 and len(sm2s) == 0:
		fc1 = sm1s[0][0]
		fc2 = sm1s[0][0]
		gamp = y[0]-sm1s[0][1]
	else:
		fc1 = default1
		fc2 = default2

	dmax = np.max(y)
	dmin = np.min(y)
	amp = dmax-dmin

	lbounds2 = [0,2500,amp/3,4,2500,amp/3,4]
	ubounds2 = [3e5,3200,2*amp,50,3200,2*amp,50]
	guess = [y[0], fc1, gamp, gwidth, fc2, gamp, gwidth]

	try:
		if (gauss):
			popt, pcov = curve_fit(fit_gaussian, x, y, p0=guess, bounds=(lbounds2,ubounds2))
#            popt, pcov = curve_fit(fit_gaussian, x, y, p0=guess)
		else:
			popt, pcov = curve_fit(fit_lorentzian, x, y, p0=guess, bounds=(lbounds2,ubounds2))
#        popt, pcov = curve_fit(fit_lorentzian, x, y, p0=guess)
#        popt, pcov = curve_fit(func, x, y, p0=guess)
	except:
		popt = [0, 0, 0, 10, 1e3, 0, 10]
		pcov = np.zeros((7,7))
		print('fit fail on file '+str(filenum))

	if (gauss):
		fit = fit_gaussian(x, *popt)
		fitg = fit_gaussian(x, *guess)
	else:
		fit = fit_lorentzian(x, *popt)
		fitg = fit_lorentzian(x, *guess)

	return popt, pcov, fit, fitg, np.transpose(mintab)