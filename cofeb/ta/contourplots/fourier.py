# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 07:53:18 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.patches as patches
import numpy as np
import load_scan as ls
import fourier_image as fi

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}
matplotlib.rc('font', **font)

scandata = ls.load_contour('/Users/alec/UCSB/scan_data/809/000809.scan')
scan_size = ls.get_scan_size('/Users/alec/UCSB/scan_data/809/000809.info')

res = len(scandata[0])

fscandata = fi.fourier_image(scandata)
fscan_size = res/scan_size

plt.close('all')
    
plt.figure(1,[5,5])
ax1 = plt.axes()

plt.imshow(np.abs(fscandata), cmap='bone', interpolation='nearest', extent=[-fscan_size/2,fscan_size/2,-fscan_size/2,fscan_size/2],
           clim=[0,1.5e5])
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)
ax1.add_patch(patches.Wedge((0, 0),6.65,0,360,fill=False,edgecolor="#ED1035" ,linewidth=2))
plt.tight_layout()
fig = plt.gcf()
fig.canvas.manager.window.raise_()

#    plt.figure(2,[6,6])    
#    plt.hist(scandata.ravel(), bins=256, range=(-1300.0, 0.0), fc='k', ec='k')

#pylab.savefig('/Users/alec/UCSB/scan_data/images/f809.pdf',bbox_inches='tight')
