# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:19:24 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import load_scan as ls

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

scanfolder = '/Users/alec/UCSB/scan_data/'

def plotcontour(scan_num,size,save=False,savefolder='/Users/alec/UCSB/scan_images/contours/',crop=[0,1,0,1],cb=True,filetype='png'):
    
    scan_num = str(scan_num)
    datapath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.scan'
    infopath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.info'
    savepath = savefolder+scan_num+'.'+filetype
    scandata = ls.load_contour(datapath)
    scan_size = ls.get_scan_size(infopath)
    lsd = len(scandata[0])
    scandata = scandata[lsd*crop[0]:lsd*crop[1],lsd*crop[2]:lsd*crop[3]]
    scan_size = scan_size*(crop[1]-crop[0])
    
    mn = np.min(scandata)
    mx = np.max(scandata)
    scandata = (scandata-mn)/(mx-mn)
    
    plt.close('all')
    
    plt.figure(1,[size,size])
    ax1 = plt.axes()
    plt.imshow(scandata, cmap='bone', interpolation='nearest')
    if cb:
        cbar = plt.colorbar(fraction=0.045, pad=0.06)
        cbar.set_ticks([0,1])
        cbar.set_ticklabels([0,1])
#    plt.axis('off')
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    plt.tight_layout()
    fig = plt.gcf()
    fig.canvas.manager.window.raise_()
    
#    plt.figure(2,[6,6])    
#    plt.hist(scandata.ravel(), bins=256, range=(-1300.0, 0.0), fc='k', ec='k')
    
    if save:
        pylab.savefig(savepath,bbox_inches='tight')

