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

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

#material parameters
pi = np.pi
surfmag = 8.336e5
sfieldpre = (1e-3)*surfmag
sres = 0.01
zfield = 7.0
zfields = zfield/sfieldpre
theta = 54.7*pi/180
phi = 180*(pi/180)

#file constants
hnum = 1
rnum = 3
dres = 50
dsize = 0.6*5
filenum = 1739

basepath = '/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/'

height = 80
blochfms = glob.glob(basepath+'bloch*_'+str(height)+'_30._*.dat')
nleftfms = glob.glob(basepath+'neelleft*_'+str(height)+'_30.*.dat')
nrightfms = glob.glob(basepath+'neelright*_'+str(height)+'_30.*.dat')

#print(blochfms)

blochcuts = []
nleftcuts = []
nrightcuts = []

#
#linecuts = np.zeros()
#
for file in blochfms:
    blochcuts.append(np.loadtxt(file))
for file in nleftfms:
    nleftcuts.append(np.loadtxt(file))
for file in nrightfms:
    nrightcuts.append(np.loadtxt(file))

blochnv = []
nleftnv = []
nrightnv = []

rad = 2

for i in range(0,hnum):
    fi = i*rnum+rad
    tlen = len(blochcuts[fi])
    blochnv.append(np.zeros((2,tlen)))
    nleftnv.append(np.zeros((2,tlen)))
    nrightnv.append(np.zeros((2,tlen)))
    simx = np.multiply(np.arange(0,sres*tlen,sres),1.1)
    blochnv[i][0] = simx
    blochnv[i][1] = np.multiply(sfieldpre,np.abs((blochcuts[fi]*np.sin(theta)*np.cos(phi))+
        (blochcuts[fi+rnum*hnum]*np.sin(theta)*np.sin(phi))+
        (np.add(blochcuts[fi+2*rnum*hnum],zfields)*np.cos(theta))))
    nleftnv[i][0] = simx
    nleftnv[i][1] = np.multiply(sfieldpre,np.abs((nleftcuts[fi]*np.sin(theta)*np.cos(phi))+
        (nleftcuts[fi+rnum*hnum]*np.sin(theta)*np.sin(phi))+
        (np.add(nleftcuts[fi+2*rnum*hnum],zfields)*np.cos(theta))))
    nrightnv[i][0] = simx
    nrightnv[i][1] = np.multiply(sfieldpre,np.abs((nrightcuts[fi]*np.sin(theta)*np.cos(phi))+
        (nrightcuts[fi+rnum*hnum]*np.sin(theta)*np.sin(phi))+
        (np.add(nrightcuts[fi+2*rnum*hnum],zfields)*np.cos(theta))))   
      
ffdata = ls.load_ff_linecut('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,7)
line = 28  
ffycut = [np.add(np.arange(0,dsize,dsize/dres),0),ffdata[0][0:dres,line],ffdata[1][0:dres,line]]
#ffycut = [np.arange(0,100),ffdata[0][:,3],ffdata[1][:,3]]           

plt.close('all')

plt.figure(1,[5,4])

imgplot = plt.imshow(ffdata[0], cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.027, pad=0.04)
plt.show()
plt.tight_layout()
#pylab.savefig('ff_1592x.pdf',bbox_inches='tight')

fig = plt.gcf()
fig.canvas.manager.window.raise_()


ploth = 0
plotr = 2
plt.figure(2,[5,4])
plt.plot(blochnv[ploth][0],blochnv[ploth][1],color='#2D7DD2',linewidth=2.0,label="Bloch")
plt.plot(nleftnv[ploth][0],nleftnv[ploth][1],color='#F97304',linewidth=2.0, label=u'left-handed Néel')
plt.plot(nrightnv[ploth][0],nrightnv[ploth][1],color='#97CC04',linewidth=2.0, label=u'right-handed Néel')
plt.errorbar(ffycut[0],ffycut[1],yerr=ffycut[2],color='#ED1035',fmt='.',label="data")
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
pylab.ylim([0,55])
plt.tight_layout()
#pylab.savefig('/Users/alec/UCSB/scan_data/images/noaxes/linecut_'+str(filenum)+'.pdf')
#
fig = plt.gcf()
fig.canvas.manager.window.raise_()

#ploth = 1
#plt.figure(2,[5,3])
#plt.fill_between(blochnv[ploth][0],blochnv[0][1],blochnv[2][1],color='#2D7DD2',alpha=0.5,linewidth=1.0)
#plt.fill_between(nleftnv[ploth][0],nleftnv[0][1],nleftnv[2][1],color='#97CC04',alpha=0.5,linewidth=1.0)
#plt.fill_between(nrightnv[ploth][0],nrightnv[0][1],nrightnv[2][1],color='#F97304',alpha=0.5,linewidth=1.0)
#pylab.savefig('linecut_1592_fillx.pdf')
#
#fig = plt.gcf()
#fig.canvas.manager.window.raise_()


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