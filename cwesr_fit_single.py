# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 07:34:44 2016

@author: alec
"""

from scipy.optimize import curve_fit
import numpy as np
import peakutils
import peakdet

gamp = 1.5e4
gwidth = 7
lbounds2 = [0,2500,1e3,4,2500,1e3,4]
ubounds2 = [3e5,3200,2e5,50,3200,2e5,50]
pdheight = 1500
maxcenshift = 20
defaultf1 = 2855
defaultf2 = 2886

def func(x, *params):
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

def cwesr_fit(x,y,filenum=0):
    [maxtab, mintab]=peakdet.peakdet(y, 7000, x)
#    indexes = peakutils.indexes(-y, thres=0.4, min_dist=4)
    sm1=[]
    sm2=[]
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
    elif len(sm1s) == 0 and len(sm2s) >= 1:
        fc1 = sm2s[0][0] 
        fc2 = sm2s[0][0]
    elif len(sm1s) >= 1 and len(sm2s) == 0:
        fc1 = sm1s[0][0]
        fc2 = sm1s[0][0]
    else:
        fc1 = defaultf1
        fc2 = defaultf2
    
    #if len(sm1s) >= 1 and len(sm2s) >= 1:
    #    guess = [edata[1,1], fc1, gamp, gwidth, fc2, gamp, gwidth]
    #elif len(sm1s) >= 0 and len(sm2s) >= 1:
    #    guess = [edata[1,1], fc1, gamp, gwidth]
    #elif len(sm1s) >= 1 and len(sm2s) >= 0:
    #    guess = [edata[1,1], fc1, gamp, gwidth]
    #else:

    guess = [y[0], fc1, gamp, gwidth, fc2, gamp, gwidth]
#    guess = [y[1], fc1, gamp, gwidth]
    
    try:
    #    if len(sm) >= 4:
    #        popt, pcov = curve_fit(func, x, y, p0=guess, bounds=(lbounds4,ubounds4))
    #    elif len(sm) == 3:
    #        popt, pcov = curve_fit(func, x, y, p0=guess, bounds=(lbounds3,ubounds3))  
    #    if len(sm1s) >= 1 and len(sm2s) >= 1:
    #        popt, pcov = curve_fit(func, x, y, p0=guess, bounds=(lbounds2,ubounds2))
    #    elif len(sm1s) == 0 and len(sm2s) >= 1:
    #        popt, pcov = curve_fit(func, x, y, p0=guess, bounds=(lbounds1,ubounds1))
    #    elif len(sm1s) >= 1 and len(sm2s) == 0:
    #        popt, pcov = curve_fit(func, x, y, p0=guess, bounds=(lbounds1,ubounds1))
    #    else:
        popt, pcov = curve_fit(func, x, y, p0=guess, bounds=(lbounds2,ubounds2))
#        popt, pcov = curve_fit(func, x, y, p0=guess)
    except:
        popt = [0, 0, 0, 10, 1e3, 0, 10]
        pcov = np.zeros((7,7))
        print('fit fail on file '+str(filenum))
        
    fit = func(x, *popt)
    fitg = func(x, *guess)
   
    return popt, pcov, fit, fitg, np.transpose(mintab)