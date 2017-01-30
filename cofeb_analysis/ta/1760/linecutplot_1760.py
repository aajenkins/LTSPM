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
import load_scan as ls
import format_plot as fp

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)

#material parameters
pi = np.pi
surfmag = 6.496e5
sfieldpre = (1e-3)*surfmag
sres = 0.01
zfield = 9.5
zfields = zfield/sfieldpre
theta = 55.60*pi/180
phi = 0*(pi/180)

#file constants
hnum = 1
rnum = 1
dres = 50
dsize = 0.6*5
filenum = 1760
rad = 400

basepath = '/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/'

heights = [68,72,76]
blochfms = np.array([])
nleftfms = np.array([])
nrightfms = np.array([])
for i in range(0,3):
    blochfms=np.append(blochfms,glob.glob(basepath+'bloch*'+str(heights[i])+'_30._'+str(rad)+'*.dat'))
    nleftfms=np.append(nleftfms,glob.glob(basepath+'neelleft*'+str(heights[i])+'_30._'+str(rad)+'*.dat'))
    nrightfms=np.append(nrightfms,glob.glob(basepath+'neelright*'+str(heights[i])+'_30._'+str(rad)+'*.dat'))


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
hnum = 3

for i in range(0,3):
    fi = 3*i
    tlen = len(blochcuts[fi])
    blochnv.append(np.zeros((2,tlen)))
    nleftnv.append(np.zeros((2,tlen)))
    nrightnv.append(np.zeros((2,tlen)))
    simx = np.add(np.multiply(np.arange(0,sres*tlen,sres),0.90),0.04)
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
      
ffdata = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,maxfgrad=15)
yconstant = 29
xconstant = 29


cutcrop = [0,50]
#
# ffycut = [np.arange(0,(cutcrop[1]-cutcrop[0])*dsize/dres,dsize/dres),ffdata[0][yline,cutcrop[0]:cutcrop[1]],ffdata[1][yline,cutcrop[0]:cutcrop[1]]]
ffxcut = [np.arange(0,(cutcrop[1]-cutcrop[0])*dsize/dres,dsize/dres),ffdata[0][yconstant,cutcrop[0]:cutcrop[1]],ffdata[1][yconstant,cutcrop[0]:cutcrop[1]]]

ffcut = ffxcut

x0, x1, y0, y1 = cutcrop[0], cutcrop[1], yconstant, yconstant


# ---------------------- PLOTS -------------------------------------------------------
# ------------------------------------------------------------------------------------

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(ffdata[0], cmap='bone', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
plt.plot([x0, x1], [y0, y1], 'r-')
plt.axis('image')
fp.format_plot(plt, 450, 450, 450, 50)

ploth = 1
plotr = 0
fig1, ax1 = plt.subplots()
plt.plot(blochnv[ploth][0],blochnv[ploth][1],color='#2D7DD2',linewidth=2.0,label="Bloch")
plt.plot(nleftnv[ploth][0],nleftnv[ploth][1],color='#F97304',linewidth=2.0, label=u'left-handed Néel')
plt.plot(nrightnv[ploth][0],nrightnv[ploth][1],color='#97CC04',linewidth=2.0, label=u'right-handed Néel')
plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="data")
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
pylab.ylim([0,40])
pylab.xlim([0,2.5])
fp.format_plot(plt, 450, 450, 0, 50)
# #pylab.savefig('/Users/alec/UCSB/scan_data/images/noaxes/linecut_'+str(filenum)+'.pdf')
# #
# #fig = plt.gcf()
# #fig.canvas.manager.window.raise_()
#
# #ploth = 1
# plt.fill_between(blochnv[ploth][0],blochnv[0][1],blochnv[2][1],color='#2D7DD2',alpha=0.5,linewidth=1.0)
# #plt.fill_between(nleftnv[ploth][0],nleftnv[0][1],nleftnv[2][1],color='#F97304',alpha=0.5,linewidth=1.0)
# plt.fill_between(nrightnv[ploth][0],nrightnv[0][1],nrightnv[2][1],color='#97CC04',alpha=0.5,linewidth=1.0)
# plt.xlabel(r'$x \quad (\mu m)$')
# plt.ylabel(r'$B_{NV} \quad (G)$')
# plt.tight_layout()
# #pylab.savefig('/Users/alec/UCSB/scan_data/images/linecut_lrgunc_'+str(filenum)+'.pdf')
#
# fig = plt.gcf()
# fig.canvas.manager.window.raise_()
#
#
# #
# #plt.figure(1,[5,4])
# #
# #imgplot = plt.imshow(ffimg, cmap='gray', interpolation='nearest',extent=[-750,750,-750,750])
# #plt.xlabel('x (nm)')
# #plt.ylabel('y (nm)')
# #plt.colorbar(fraction=0.046, pad=0.04)
# #plt.show()
# #plt.tight_layout()
# #pylab.savefig('ff1.pdf')

plt.show()