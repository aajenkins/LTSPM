# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:17:32 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import glob
import re
import load_scan as ls

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)

#material parameters
pi = np.pi
surfmag = 7.393e5
sfieldpre = (1e-3)*surfmag
sres = 0.01
zfield = 0
zfields = zfield/sfieldpre
theta = 55.76*pi/180
phi = 180*(pi/180)

#file constants
hnum = 1
rnum = 1
dres = 60
vsize = 0.8
dsize = vsize*5
filenum = 1847
rad = 400

basepath = '/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/'

heights = [70,80,90]
blochfms = np.array([])
nleftfms = np.array([])
nrightfms = np.array([])
for i in range(0,3):
    blochfms=np.append(blochfms,glob.glob(basepath+'bloch*_'+str(heights[i])+'_30._'+str(rad)+'*.dat'))
    nleftfms=np.append(nleftfms,glob.glob(basepath+'neelleft*_'+str(heights[i])+'_30._'+str(rad)+'*.dat'))
    nrightfms=np.append(nrightfms,glob.glob(basepath+'neelright*_'+str(heights[i])+'_30._'+str(rad)+'*.dat'))
print(blochfms)

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
    simx = np.add(np.multiply(np.arange(0,sres*tlen,sres),0.75),0.15)
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
yconstant = 26
xconstant = 27

cutcrop = [10,48]

ffxcut = [np.arange(0,(cutcrop[1]-cutcrop[0])*dsize/dres,dsize/dres),ffdata[0][yconstant,cutcrop[0]:cutcrop[1]],ffdata[1][yconstant,cutcrop[0]:cutcrop[1]]]
ffycut = [np.arange(0,(cutcrop[1]-cutcrop[0])*dsize/dres,dsize/dres),ffdata[0][cutcrop[0]:cutcrop[1],xconstant],ffdata[1][cutcrop[0]:cutcrop[1],xconstant]]

ffcut = ffxcut
x0, x1, y0, y1 = cutcrop[0], cutcrop[1], xconstant, xconstant

#ffcut = ffycut
#x0, x1, y0, y1 = yconstant, yconstant, 1, dres-1

plt.close('all')

plt.figure(1,[5,5])

imgplot = plt.imshow(ffdata[0], cmap='bone', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.plot([x0, x1], [y0, y1], 'r-')
plt.axis('image')
plt.show()
plt.tight_layout()
#pylab.savefig('/Users/alec/UCSB/scan_data/images/ff_'+str(filenum)+'.pdf',bbox_inches='tight')

fig = plt.gcf()
fig.canvas.manager.window.raise_()


ploth = 1
plotr = 0
plt.figure(2,[8,6])
plt.plot(blochnv[ploth][0],blochnv[ploth][1],color='#2D7DD2',linewidth=2.0,label="Bloch")
#plt.plot(nleftnv[ploth][0],nleftnv[ploth][1],color='#F97304',linewidth=2.0, label=u'left-handed Néel')
plt.plot(nleftnv[ploth][0],nleftnv[ploth][1],color='#97CC04',linewidth=2.0, label=u'left-handed Néel')
plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="data")
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
pylab.ylim([0,50])
#pylab.xlim([0,2.5])
plt.tight_layout()
#pylab.savefig('/Users/alec/UCSB/scan_data/images/noaxes/linecut_'+str(filenum)+'.pdf')
#
#fig = plt.gcf()
#fig.canvas.manager.window.raise_()

#ploth = 1
plt.fill_between(blochnv[ploth][0],blochnv[0][1],blochnv[2][1],color='#2D7DD2',alpha=0.5,linewidth=1.0)
#plt.fill_between(nleftnv[ploth][0],nleftnv[0][1],nleftnv[2][1],color='#F97304',alpha=0.5,linewidth=1.0)
plt.fill_between(nleftnv[ploth][0],nleftnv[0][1],nleftnv[2][1],color='#97CC04',alpha=0.5,linewidth=1.0)
plt.xlabel(r'$x \quad (\mu m)$')
plt.ylabel(r'$B_{NV} \quad (G)$')
plt.tight_layout()
#pylab.savefig('/Users/alec/UCSB/scan_data/images/linecut_lrgunc_'+str(filenum)+'.pdf')

fig = plt.gcf()
fig.canvas.manager.window.raise_()


#
#plt.figure(1,[5,4])
#
#imgplot = plt.imshow(ffimg, cmap='gray', interpolation='nearest',extent=[-750,750,-750,750])
#plt.xlabel('x (nm)')
#plt.ylabel('y (nm)')
#plt.colorbar(fraction=0.046, pad=0.04)
#plt.show()
#plt.tight_layout()
#pylab.savefig('ff1.pdf')