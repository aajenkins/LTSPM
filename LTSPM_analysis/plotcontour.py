# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:19:24 2016

@author: alec
"""

# import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import load_scan as lscan
import plotting.format_plots_tkagg as fp
from PIL import Image

scan_num = 824
size = 10e-6
save=False
savefolder='/Users/alec/UCSB/scan_images/contours/'
filetype='png'
scanfolder = '/Users/alec/UCSB/scan_data/'

# def plotcontour(scan_num, size, save=False, savefolder='/Users/alec/UCSB/scan_images/contours/', crop=[0,1,0,1], filetype='png'):

scan_num = str(scan_num)
datapath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.scan'
infopath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.info'
savepath = savefolder+scan_num+'.'+filetype
scandata = lscan.load_contour(datapath)[1]
scan_size = lscan.get_scan_size(infopath)
mn = np.min(scandata)
mx = np.max(scandata)
# scandata = (scandata-mn)/(mx-mn)

# scandata = (1e-3)*scandata

plt.close('all')
my_dpi = 96
size_inches = size/my_dpi
if save:
    fig = plt.figure(frameon=False)
    fig.set_size_inches(size_inches, size_inches)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    im = plt.imshow(scandata, cmap='bone', interpolation="nearest")
    # cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    # cbar.set_ticks([0,1])
    # cbar.set_ticklabels([0,1])
    # pcol.axes.set_xlim(0, len(scandata[0]))
    # pcol.axes.set_ylim(0, len(scandata))
    plt.gca().invert_yaxis()
    ax.set_axis_off()
    fig.savefig(savepath, format=filetype, dpi=my_dpi)
else:
    fig, ax = plt.subplots()
    im = plt.imshow(scandata, cmap='bone', interpolation='nearest')
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    # cbar.set_ticks([0,1])
    # cbar.set_ticklabels([0,1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    cbar.set_label('NV PL difference: RF on - RF off (kcnts)')
    plt.subplots_adjust(right=0.8)

plt.show()

# if __name__ == "__main__":
#     import sys
#     if (len(sys.argv) == 3):
#         plot_stray_field(int(sys.argv[1]), eval(sys.argv[2]))
#     else:
#         print('enter scan number and scan size')
