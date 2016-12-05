# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 11:00:39 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import load_scan as ls
import matplotlib.patches as patches

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

scanfolder = '/Users/alec/UCSB/scan_data/'
savefolder = '/Users/alec/UCSB/scan_data/images/noaxes/'

scan_num = 135

scan_num = str(scan_num)
datapath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.scan'
infopath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.info'
savepath = savefolder+scan_num+'.pdf'
scandata = ls.load_contour(datapath)
scan_size = ls.get_scan_size(infopath)

scandata = scandata/(1e3)

plt.close('all')

fp = plt.figure(1,[4,4])
ax1 = plt.axes()

plt.imshow(scandata, cmap='bone', interpolation='nearest', extent=[-scan_size/2,scan_size/2,-scan_size/2,scan_size/2],clim=[44,62])
cbar = plt.colorbar(fraction=0.045, pad=0.06)
cbar.set_ticks([44,62])
cbar.set_ticklabels([44,62])
#plt.axis('off')
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)
plt.tight_layout()
fig = plt.gcf()
fig.canvas.manager.window.raise_()

#plt.figure(2,[6,6])    
#plt.hist(scandata.ravel(), bins=256, range=(0, 7e4), fc='k', ec='k')

pylab.savefig(savepath,bbox_inches='tight')