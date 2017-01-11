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
from scipy import misc
import matplotlib.pylab as pylab
import load_scan as ls
import fourier_image as fi

pi = np.pi

scannum = 1760
#scanbacknum = 1740
xres = 50
yres = 50
zfield = 9.5
scanL = 0.6*5000
height = 66


data = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,10)
#misc.imsave('/Users/alec/UCSB/scan_data/images/ff1760.png', data[0])
#data_back = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(scanbacknum)+'-esrdata/fitdata.txt',xres,yres,7)
ffmask = ndimage.imread('/Users/alec/UCSB/scan_data/images/ff1760mask_blur.png',flatten=True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

#---------------- FIT FUNCTIONS ----------------------------------
#-----------------------------------------------------------------

def fitarctan(x, *params):
    y = np.zeros_like(x)
    c = params[0]
    a = params[1]
    x0 = params[2]
    wid = params[3]
    
    y = c+a*np.arctan((x-x0)/wid)
    return y

#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

phi = 0
theta = 55.76*np.pi/180

datas = np.multiply(ffmask,data[0])
datas0 = np.add(datas,-np.cos(theta)*zfield)
#datasback = data_back[0]
#dbackshift = np.concatenate((np.full((3,49),zfield*np.cos(theta)),(datasback[:-3,1:])),axis=0)
#dbackshift = np.concatenate((dbackshift,np.full((50,1),zfield*np.cos(theta))),axis=1)
#datadiff = datas-dbackshift

#datas = np.multiply(1,data[0])
wimdata = fi.window_image(datas0)
fdata = fi.fourier_image(datas0)

dlen = len(fdata)
hlen = int(np.floor(dlen/2))

hzf = np.zeros_like(fdata)
k = 0

for j in range(0,dlen):
    for i in range(0,dlen):
        rij = np.sqrt((i-hlen)**2+(j-hlen)**2)
        k = 2*pi*rij/(scanL)
        if (i==hlen and j==hlen):
            hzf[j,i] = fdata[j,i]/np.sin(theta)
        else:
            hzf[j,i] = fdata[j,i]/(np.cos(theta)*(1-
            1j*((i-hlen)/rij)*np.tan(theta)))

hz = np.real(fi.ifourier_image(hzf))

#---------------- LINECUTS ---------------------------------------
#-----------------------------------------------------------------

x0, y0 = 23.5, 29
phinum = 8
lcnum = 12
lclen = 12
hzphi = np.zeros((phinum,lcnum))
for i in range(0,phinum):    
    phi = i*pi/4
    x1, y1 = x0-lclen*np.cos(phi), y0-lclen*np.sin(phi)
    x, y = np.linspace(x0, x1, lcnum), np.linspace(y0, y1, lcnum)
    hzphi[i] = ndimage.map_coordinates(np.transpose(hz), np.vstack((x,y)))
    


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
fig.setGeometry(50,50,dx, dy)
fig.raise_()


plt.figure(3,[4.5,4.5])
im = plt.imshow(hz, cmap='jet', interpolation='nearest')
plt.colorbar(im, fraction=0.046, pad=0.04, use_gridspec=True)

for i in range(0,phinum):
    phi = i*pi/4
    x1, y1 = x0-lclen*np.cos(phi), y0-lclen*np.sin(phi)
    plt.plot([x0, x1], [y0, y1], 'ro-')
plt.axis('image')

plt.tight_layout()

fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(420,50,dx, dy)
fig.raise_()


fig, axes = plt.subplots(nrows=phinum, sharex=True, sharey=True, figsize=(5,9))
for i in range(0,phinum):
    axes[i].plot(hzphi[i])
    axes[i].get_yaxis().set_visible(False)

fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
#fig.tight_layout()

fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(800,50,dx, dy)
fig.raise_()


plt.figure(6,[4.5,4.5])
plt.imshow(hz, cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(50,450,dx, dy)
fig.raise_()
#pylab.savefig('/Users/alec/UCSB/scan_data/images/reconstruction/ff_recon_bz_'+str(filenum)+'.pdf')
