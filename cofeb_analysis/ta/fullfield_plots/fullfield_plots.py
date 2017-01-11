# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 23:26:30 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
from scipy import ndimage
import load_scan as ls

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

scannum = 1739
bscannum = 1740
hnum = 1
rnum = 3
dres = 50
dsize = 0.6*5


xres = 50
yres = 50   
#ffdata = ls.load_rf_track('/Users/alec/UCSB/scan_data/'+str(scannum)+'/'+str(scannum).zfill(6)+'.scan',xres,yres)
ffdata = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,7)
ffback = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(bscannum)+'-esrdata/fitdata.txt',xres,yres,7)
ffdatac = ffdata[0][1:50,0:48]
sub = np.zeros((5, 5))
#for j in range(0,5):
#    for i in range(0,5):
#        ffbackc = ffback[0][j:j+45,i+1:i+46]
#        ffdatasub = np.add(ffdatac,np.multiply(ffbackc,-1))
#        ffdatacmask = ffdatasub
#        ffdatacmask[0:14,0:11] = np.zeros_like(ffdatacmask[0:14,0:11])
#        ffdatacmask[10:35,8:33]= np.zeros_like(ffdatacmask[10:35,8:33])
#        sub[i,j] = np.abs(np.sum(ffdatacmask))

ffbackc = ffback[0][0:49,2:50]
ffdatasub = np.add(ffdatac,np.multiply(ffbackc,-1))

#ffmask = ndimage.imread('/Users/alec/UCSB/scan_data/images/1739ffmask.png')
plt.close('all')

plt.figure(1,[4,4])
ax1 = plt.axes()

imgplot = plt.imshow(ffdata[0], cmap='bone', interpolation='None')
ax1.axis('off')
#plt.colorbar(fraction=0.045, pad=0.04)
#ax1.get_xaxis().set_visible(False)
#ax1.get_yaxis().set_visible(False)
#plt.tight_layout()
pylab.savefig('/Users/alec/UCSB/scan_data/images/'+str(scannum)+'ff.png',bbox_inches='tight')

fig = plt.gcf()
fig.canvas.manager.window.raise_()

#plt.figure(2,[4,4])
#ax1 = plt.axes()
#imgplot = plt.imshow(ffdatac, cmap='bone', interpolation='None')
##plt.imshow(ffbackc, cmap='bone', alpha=0.9, interpolation='nearest')
#plt.colorbar(fraction=0.045, pad=0.04)
#ax1.get_xaxis().set_visible(False)
#ax1.get_yaxis().set_visible(False)
#plt.tight_layout()
##pylab.savefig('/Users/alec/UCSB/scan_data/images/'+str(scannum)+'ff.pdf',bbox_inches='tight')
#
#fig = plt.gcf()
#fig.canvas.manager.window.raise_()

#line = 28
#ffycut = [np.add(np.arange(0,dsize,dsize/dres),0.65),ffdata[0][line,0:dres],ffdata[1][line,0:dres]]    
#
#ploth = 0
#plotr = 2
#plt.figure(2,[5,4])
##plt.plot(blochnv[ploth][0],blochnv[ploth][1],color='#2D7DD2',linewidth=2.0,label="Bloch")
##plt.plot(nleftnv[ploth][0],nleftnv[ploth][1],color='#F97304',linewidth=2.0, label=u'left-handed Néel')
##plt.plot(nrightnv[ploth][0],nrightnv[ploth][1],color='#97CC04',linewidth=2.0, label=u'right-handed Néel')
#plt.errorbar(ffycut[0],ffycut[1],yerr=ffycut[2],color='#ED1035',fmt='.',label="data")
#plt.legend(loc=1,borderaxespad=1,prop={'size':10})
#pylab.ylim([0,55])
#plt.tight_layout()
##pylab.savefig('/Users/alec/UCSB/scan_data/images/noaxes/linecut_717.pdf')
##
#fig = plt.gcf()
#fig.canvas.manager.window.raise_()
