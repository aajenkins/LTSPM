# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:17:32 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy import misc
import numpy as np
import glob
import load_scan as ls
import format_plot as fp

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)

#material parameters
pi = np.pi
surfmag = 1.468e6
sfieldpre = (1e-3)*surfmag
sres = 0.01
zfield = 0
zfields = zfield/sfieldpre
theta = 55.06*pi/180
phi = 180*(pi/180)

#file constants
hnum = 1
rnum = 1
dres = 50
vsize = 0.5
dsize = vsize*5
filenum = 1903
rad = 400

basepath = '/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/'

heights = [134,140,146]
blochfms = np.array([])
nleftfms = np.array([])
nrightfms = np.array([])
for i in range(0,3):
    blochfms=np.append(blochfms,glob.glob(basepath+'bloch*_'+str(heights[i])+'_30._'+str(rad)+'*.dat'))
    nleftfms=np.append(nleftfms,glob.glob(basepath+'neelleft*_'+str(heights[i])+'_30._'+str(rad)+'*.dat'))
    nrightfms=np.append(nrightfms,glob.glob(basepath+'neelright*_'+str(heights[i])+'_30._'+str(rad)+'*.dat'))

blochcuts = []
nleftcuts = []
nrightcuts = []

for file in blochfms:
    blochcuts.append(np.loadtxt(file))
for file in nleftfms:
    nleftcuts.append(np.loadtxt(file))
for file in nrightfms:
    nrightcuts.append(np.loadtxt(file))

blochnv = []
nleftnv = []
nrightnv = []

rad = 0
hnum = len(heights)

for i in range(0,hnum):
    fi = 3*i
    tlen = len(blochcuts[fi])
    blochnv.append(np.zeros((2,tlen)))
    nleftnv.append(np.zeros((2,tlen)))
    nrightnv.append(np.zeros((2,tlen)))
    simx = np.add(np.multiply(np.arange(0,sres*tlen,sres),0.84),-0.15)
    blochnv[i][0] = simx
    blochnv[i][1] = np.multiply(sfieldpre,np.abs((blochcuts[fi]*np.sin(theta)*np.cos(phi))+
        (blochcuts[fi+1]*np.sin(theta)*np.sin(phi))+
        (np.add(blochcuts[fi+2],zfields)*np.cos(theta))))
    nleftnv[i][0] = simx
    nleftnv[i][1] = np.multiply(sfieldpre,np.abs((nleftcuts[fi]*np.sin(theta)*np.cos(phi))+
        (nleftcuts[fi+1]*np.sin(theta)*np.sin(phi))+
        (np.add(nleftcuts[fi+2],zfields)*np.cos(theta))))
    nrightnv[i][0] = simx
    nrightnv[i][1] = np.multiply(sfieldpre,np.abs((nrightcuts[fi]*np.sin(theta)*np.cos(phi))+
        (nrightcuts[fi+1]*np.sin(theta)*np.sin(phi))+
        (np.add(nrightcuts[fi+2],zfields)*np.cos(theta))))

ffdata = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,maxfgrad=20)
# misc.imsave('/Users/alec/UCSB/scan_images/full-field/ff'+str(filenum)+'.png', ffdata[0])

yconstant = 28
xconstant = 21

cutcrop = [0,50]

ffxcut = [np.arange(0,(cutcrop[1]-cutcrop[0])*dsize/dres,dsize/dres),ffdata[0][yconstant,cutcrop[0]:cutcrop[1]],ffdata[1][yconstant,cutcrop[0]:cutcrop[1]]]
ffycut = [np.arange(0,(cutcrop[1]-cutcrop[0])*dsize/dres,dsize/dres),ffdata[0][cutcrop[0]:cutcrop[1],xconstant],ffdata[1][cutcrop[0]:cutcrop[1],xconstant]]

ffcut = ffxcut
x0, x1, y0, y1 = cutcrop[0], cutcrop[1], yconstant, yconstant

# ffcut = ffycut
# x0, x1, y0, y1 = xconstant, xconstant, 1, dres-1

x0, x1, y0, y1 = np.multiply(dsize/dres,[x0, x1, y0, y1])

fig1, ax1 = plt.subplots()
im1 = plt.plot(blochcuts[2])
fp.format_plot(plt, 350, 350, 50, 450)


fig1, ax1 = plt.subplots()
im1 = plt.imshow(ffdata[0], cmap='bone', interpolation='nearest',extent=[0,dsize,dsize,0])
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
plt.plot([x0, x1], [y0, y1], 'r-')
plt.axis('image')
fp.format_plot(plt, 350, 350, 50, 50)
# pylab.savefig('/Users/alec/UCSB/scan_images/ff_'+str(filenum)+'.png')


ploth = 1
plotr = 0
fig2, ax2 = plt.subplots()
plt.plot(blochnv[ploth][0],blochnv[ploth][1],color='#2D7DD2',linewidth=2.0,label="Bloch")
plt.plot(nrightnv[ploth][0],nrightnv[ploth][1],color='#F97304',linewidth=2.0, label=u'right-handed Néel')
plt.plot(nleftnv[ploth][0],nleftnv[ploth][1],color='#97CC04',linewidth=2.0, label=u'left-handed Néel')
plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="data")
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
pylab.ylim([0,50])
plt.fill_between(blochnv[ploth][0],blochnv[0][1],blochnv[2][1],color='#2D7DD2',alpha=0.5,linewidth=1.0)
plt.fill_between(nrightnv[ploth][0],nrightnv[0][1],nrightnv[2][1],color='#F97304',alpha=0.5,linewidth=1.0)
plt.fill_between(nleftnv[ploth][0],nleftnv[0][1],nleftnv[2][1],color='#97CC04',alpha=0.5,linewidth=1.0)
plt.xlim([0,2.4])
plt.xlabel(r'$x \quad (\mu m)$')
plt.ylabel(r'$B_{NV} \quad (G)$')
fp.format_plot(plt, 500, 450, 450, 50)
pylab.savefig('/Users/alec/UCSB/scan_images/linecut_'+str(filenum)+'.png')

plt.show()