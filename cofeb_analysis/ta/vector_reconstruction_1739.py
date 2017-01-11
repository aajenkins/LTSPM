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
import matplotlib.pylab as pylab
import load_scan as ls
import fourier_image as fi

pi = np.pi

scannum = 1739
scanbacknum = 1740
xres = 50
yres = 50
zfield = 9.0
scanL = 0.6*5000
height = 66


data = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,10)
data_back = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(scanbacknum)+'-esrdata/fitdata.txt',xres,yres,7)
ffmask = ndimage.imread('/Users/alec/UCSB/scan_data/images/1739ffmask.png',flatten=True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

phi = 0
theta = 55.76*np.pi/180

datas = np.multiply(ffmask,data[0])
datas0 = np.add(datas,-np.cos(theta)*zfield)
datasback = data_back[0]
dbackshift = np.concatenate((np.full((3,49),zfield*np.cos(theta)),(datasback[:-3,1:])),axis=0)
dbackshift = np.concatenate((dbackshift,np.full((50,1),zfield*np.cos(theta))),axis=1)
datadiff = datas-dbackshift

#datas = np.multiply(1,data[0])
wimdata = fi.window_image(datas0)
fdata = fi.fourier_image(datas0)

dlen = len(fdata)
hlen = int(np.floor(dlen/2))

hzf = np.zeros_like(fdata)
mchargef = np.zeros_like(fdata)
k = 0

for j in range(0,dlen):
    for i in range(0,dlen):
        rij = np.sqrt((i-hlen)**2+(j-hlen)**2)
        k = 2*pi*rij/(scanL)
        if (i==hlen and j==hlen):
            hzf[j,i] = fdata[j,i]/np.sin(theta)
            mchargef[j,i] = 0
        else:
            hzf[j,i] = fdata[j,i]/(np.cos(theta)*(1-
            1j*((i-hlen)/rij)*np.tan(theta)))
            mchargef[j,i] = (1/(k))*np.exp(height*k)*hzf[j,i]

rdata = np.real(fi.ifourier_image(hzf))
mdata = np.real(fi.ifourier_image(mchargef))

#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

plt.figure(1,[4.5,4.5])
plt.imshow(datas, cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(50,450,dx, dy)
fig.raise_()

plt.figure(2,[4.5,4.5])
plt.imshow(datasback, cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(500,450,dx, dy)
fig.raise_()
pylab.savefig('/Users/alec/UCSB/scan_data/images/reconstruction/ff_'+str(scanbacknum)+'.pdf')

plt.figure(6,[4.5,4.5])
plt.imshow(datadiff, cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(1000,450,dx, dy)
fig.raise_()

plt.figure(3,[4.5,4.5])
plt.imshow(mdata, cmap='jet', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(50,50,dx, dy)
fig.raise_()

plt.figure(4,[4.5,4.5])
plt.plot(mdata[28,:])
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(1000,50,dx, dy)
fig.raise_()


plt.figure(5,[4.5,4.5])
plt.imshow(rdata, cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(500,50,dx, dy)
fig.raise_()
#pylab.savefig('/Users/alec/UCSB/scan_data/images/reconstruction/ff_recon_bz_'+str(filenum)+'.pdf')
