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
from scipy.fftpack import fft2, ifft2
import load_scan as ls
import fourier_image as fi

scannum = 1739
xres = 50
yres = 50

data = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,7)
ffmask = ndimage.imread('/Users/alec/UCSB/scan_data/images/1739ffmaskaa.png',flatten=True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

phi = 0
theta = -54.7*np.pi/180

datas = np.multiply(ffmask,data[0])
fdata = fi.fourier_image(datas)
wimdata = fi.window_image(datas)

dlen = len(fdata)
hlen = int(np.floor(dlen/2))

hzf = np.zeros_like(fdata)

for j in range(0,dlen):
    for i in range(0,dlen):
        if (i==hlen and j==hlen):
            hzf[j,i] = fdata[j,i]/np.sin(theta)
        else:
            hzf[j,i] = fdata[j,i]/(np.cos(theta)*(1+
            1j*((i-hlen)/np.sqrt((i-hlen)**2+(j-hlen)**2))*np.tan(theta)))
        if np.isnan(hzf[j,i]):
            print(j,i)

rdata = np.add(np.add(np.real(fi.ifourier_image(hzf)),1),-1)


#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

plt.figure(1,[4,4])
plt.imshow(np.log(np.abs(fdata)+1), cmap='gray', interpolation='nearest')
fig = plt.gcf()
fig.canvas.manager.window.raise_()
plt.figure(2,[4,4])
plt.imshow(wimdata, cmap='gray', interpolation='nearest')
fig = plt.gcf()
fig.canvas.manager.window.raise_()
plt.figure(3,[4,4])
plt.imshow(datas, cmap='gray', interpolation='nearest')
fig = plt.gcf()
fig.canvas.manager.window.raise_()
plt.figure(4,[4,4])
plt.imshow(rdata, cmap='bone', interpolation='nearest')
plt.colorbar()
fig = plt.gcf()
fig.canvas.manager.window.raise_()