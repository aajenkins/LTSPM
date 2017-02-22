# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 22:02:57 2016

@author: alec
"""

import numpy as np
from scipy.fftpack import fft2, ifft2

def fourier_image (imdata):
    wdata = window_image(imdata)
    fdata = fft2(imdata)
    cenfdata = move_quad(fdata)
    return cenfdata

def ifourier_image (fdata):
    decenfdata = move_quad(fdata)
    rdata = ifft2(decenfdata)
    return rdata

def window_image (imdata, power=(1/2)):
    wimdata = np.zeros_like(imdata)
    dlen = len(imdata)
    for j in range(0,dlen):
        for i in range(0,dlen):
            wimdata[i][j] = imdata[i][j]*((1/4)*((1-np.cos(2*np.pi*(i/dlen))))*((1-np.cos(2*np.pi*(j/dlen)))))**power
    return wimdata

def move_quad(data):
    dlen = len(data)
    hlen = int(np.floor(dlen/2))

    quad1 = data[0:hlen,0:hlen]
    quad2 = data[hlen:dlen,0:hlen]
    quad3 = data[0:hlen,hlen:dlen]
    quad4 = data[hlen:dlen,hlen:dlen]
    movedata = np.concatenate((np.concatenate((quad4,quad2),axis=1),np.concatenate((quad3,quad1),axis=1)),axis=0)

    return movedata
