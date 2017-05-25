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
import format_plot as fp

# font = {'family' : 'Arial',
#         'weight' : 'normal',
#         'size'   : 15}
# matplotlib.rc('font', **font)
# plt.style.use('ggplot')

scanfolder = '/Users/alec/UCSB/scan_data/'

def plotcontour(scan_num, size, save=False, savefolder='/Users/alec/UCSB/scan_images/contours/', crop=[0,1,0,1], filetype='png'):

    scan_num = str(scan_num)
    datapath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.scan'
    infopath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.info'
    savepath = savefolder+scan_num+'.'+filetype
    scandata = lscan.load_contour(datapath)
    scan_size = lscan.get_scan_size(infopath)
    lsd = len(scandata[0])
    scandata = scandata[int(lsd*crop[0]):int(lsd*crop[1]),int(lsd*crop[2]):int(lsd*crop[3])]
    scan_size = scan_size*(crop[1]-crop[0])
    mn = np.min(scandata)
    mx = np.max(scandata)
    scandata = (scandata-mn)/(mx-mn)

    plt.close('all')
    my_dpi = 96
    size_inches = size/my_dpi
    if save:
        fig = plt.figure(frameon=False)
        fig.set_size_inches(size_inches, size_inches)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        fig.add_axes(ax)
        im = plt.imshow(scandata, cmap='bone', interpolation="nearest")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_ticks([0,1])
        cbar.set_ticklabels([0,1])
        # pcol.axes.set_xlim(0, len(scandata[0]))
        # pcol.axes.set_ylim(0, len(scandata))
        plt.gca().invert_yaxis()
        ax.set_axis_off()
        fig.savefig(savepath, format=filetype, dpi=my_dpi)
    else:
        fig, ax = plt.subplots()
        im = plt.imshow(scandata, cmap='bone', interpolation='nearest')
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_ticks([0,1])
        cbar.set_ticklabels([0,1])
        fp.format_plot(plt, size, size, 0, 50)
