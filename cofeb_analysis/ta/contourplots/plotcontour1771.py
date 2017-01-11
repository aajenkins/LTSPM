# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:19:24 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

def load_contour (path):
    scanlist = np.loadtxt(path, delimiter='\t')[:,6]
    scan2d = np.zeros((10, 100))
#    print(scanlist)
    for j in range(0,10):
        for i in range(0,100):
            scan2d[j,i] = scanlist[j+(i*10)]
    return scan2d

scanfolder = '/Users/alec/UCSB/scan_data/'
 
scan_num = str(1771)
size = 5
datapath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.scan'
infopath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.info'
#savepath = savefolder+scan_num+'.pdf'
scandata = load_contour(datapath)
#scandata = scandata[lsd*crop[0]:lsd*crop[1],lsd*crop[2]:lsd*crop[3]]
#scan_size = scan_size*(crop[1]-crop[0])

#mn = np.min(scandata)
#mx = np.max(scandata)
#scandata = (scandata-mn)/(mx-mn)

cut = scandata[2,:]

plt.close('all')

plt.figure(1,[size,size])
ax1 = plt.axes()
plt.plot(cut)
#if cb:
#    cbar = plt.colorbar(fraction=0.045, pad=0.06)
#    cbar.set_ticks([0,1])
#    cbar.set_ticklabels([0,1])
#    plt.axis('off')
#ax1.get_xaxis().set_visible(False)
#ax1.get_yaxis().set_visible(False)
plt.tight_layout()
fig = plt.gcf()
fig.canvas.manager.window.raise_()



#    plt.figure(2,[6,6])    
#    plt.hist(scandata.ravel(), bins=256, range=(-1300.0, 0.0), fc='k', ec='k')


