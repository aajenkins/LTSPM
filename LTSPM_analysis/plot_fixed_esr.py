# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:19:24 2016

@author: alec
"""

# import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
from matplotlib_scalebar.scalebar import ScaleBar

import load_scan as lscan
import plotting.format_plots_tkagg as fp
from PIL import Image

scan_num = 1930
size = 5e-5
save=False
savefolder='/Users/alec/UCSB/scan_images/contours/'
filetype='png'
scanfolder = '/Users/alec/UCSB/scan_data/'

scan_num = str(scan_num)
datapath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.scan'
infopath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.info'
savepath = savefolder+scan_num+'.'+filetype
scandata = lscan.load_contour(datapath)[1]
scan_size = lscan.get_scan_size(infopath)
mn = np.min(scandata)
mx = np.max(scandata)

# scandata = 2*scandata

plt.close('all')

fig, ax = plt.subplots()
im = plt.imshow(scandata, cmap='bone', interpolation='nearest')
scalebar = ScaleBar(0.2)
plt.gca().add_artist(scalebar)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)
cbar.set_label('')
plt.subplots_adjust(right=0.8)

plt.show()
