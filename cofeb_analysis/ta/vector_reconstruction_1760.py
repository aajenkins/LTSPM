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
import load_scan as ls
import fourier_image as fi

pi = np.pi

scannum = 1760
xres = 50
yres = 50
height = (50/(0.6*5000))*2


data = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,10)
ffmask = ndimage.imread('/Users/alec/UCSB/scan_data/images/ff_1760mask.png',flatten=True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

phi = 0
theta = 55.76*np.pi/180

datas = np.multiply(ffmask,data[0])
datas0 = np.add(datas,-np.cos(theta)*9.5)
#datas = np.multiply(1,data[0])
wimdata = fi.window_image(datas0)
fdata = fi.fourier_image(wimdata)

dlen = len(fdata)
hlen = int(np.floor(dlen/2))

hzf = np.zeros_like(fdata)
mchargef = np.zeros_like(fdata)
k = 0

for j in range(0,dlen):
    for i in range(0,dlen):
        k = np.sqrt((i-hlen)**2+(j-hlen)**2)
        if (i==hlen and j==hlen):
            hzf[j,i] = fdata[j,i]/np.sin(theta)
            mchargef[j,i] = 0
        else:
            hzf[j,i] = fdata[j,i]/(np.cos(theta)*(1-
            1j*((i-hlen)/k)*np.tan(theta)))
            mchargef[j,i] = (1/(2*pi*k))*np.exp(2*pi*height*k)*hzf[j,i]

rdata = np.real(fi.ifourier_image(hzf))
mdata = np.real(fi.ifourier_image(mchargef))

#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

plt.figure(1,[5,5])
plt.imshow(datas, cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(50,400,dx, dy)
fig.raise_()

plt.figure(2,[5,5])
plt.imshow(datas0, cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(500,400,dx, dy)
fig.raise_()

plt.figure(3,[5,5])
plt.imshow(mdata, cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(50,50,dx, dy)
fig.raise_()

plt.figure(4,[5,5])
plt.plot(mdata[28,:])
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(900,200,dx, dy)
fig.raise_()


plt.figure(5,[5,5])
plt.imshow(rdata, cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(500,50,dx, dy)
fig.raise_()
#pylab.savefig('/Users/alec/UCSB/scan_data/images/reconstruction/ff_recon_bz_'+str(filenum)+'.pdf')
