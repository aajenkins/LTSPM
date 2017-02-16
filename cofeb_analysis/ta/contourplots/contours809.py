# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 11:00:39 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
from scipy import misc
from scipy import signal
from skimage import morphology
from skimage import exposure
import load_scan as lscan
import matplotlib.patches as patches
import format_plot as fp
from skimage import measure
from skimage import morphology
from skimage import feature

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

scanfolder = '/Users/alec/UCSB/scan_data/'
savefolder = '/Users/alec/UCSB/scan_data/images/noaxes/'

scan_num = 809

scan_num = str(scan_num)
datapath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.scan'
infopath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.info'
savepath = savefolder+scan_num+'.pdf'
scandata = lscan.load_contour(datapath)
scan_size = lscan.get_scan_size(infopath)

dmin = np.min(scandata)
dmax = np.max(scandata)
contrast = 1.5
slen = len(scandata)

scandata = np.multiply(np.add(scandata,-dmin),1/(dmax-dmin))
binaryScanData = scandata.copy()
# scandata = exposure.adjust_sigmoid(scandata)
binary_cutoff = 0.5
for i in range(slen):
    for j in range(slen):
        if (binaryScanData[j,i] < binary_cutoff):
            binaryScanData[j,i] = 0
        else:
            binaryScanData[j,i] = 1

scandataopen = morphology.binary_opening(binaryScanData, morphology.disk(1))
skeleton = morphology.skeletonize(binaryScanData)

# binarydata = scandata.copy()
#
# for j in range(0, slen):
#     for i in range(0, slen):
#         if (binarydata[j,i] >= 0.55):
#             binarydata[j,i] = 0
#         else:
#             binarydata[j,i] = 1
#
# selem1=morphology.disk(1)
# selem3=morphology.disk(3)
#
# # scandata_close = morphology.closing(scandata, selem)
# binarydata_close =  morphology.binary_closing(binarydata, selem1)
# skeleton = morphology.skeletonize(binarydata_close)
# skeleton_close = morphology.binary_closing(skeleton,selem3)
# doubleskeleton = morphology.skeletonize(skeleton_close)
# contours = measure.find_contours(scandata, 0.5)

# for n, contour in enumerate(contours):
#     if (len(contour)<20):

misc.imsave('/Users/alec/UCSB/scan_images/contours/'+str(scan_num)+'.png', scandata)

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(binaryScanData, cmap='Greys', interpolation='nearest')
fp.format_plot(plt, 400, 400, 50, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(skeleton, cmap='Greys', interpolation='nearest')
fp.format_plot(plt, 400, 400, 450, 50)

# fig1, ax1 = plt.subplots()
# im1 = plt.imshow(skeleton, cmap='Greys', interpolation='nearest')
# fp.format_plot(plt, 400, 400, 450, 50)
#
# fig1, ax1 = plt.subplots()
# im1 = plt.imshow(binarydata, cmap='Greys', interpolation='nearest')
# fp.format_plot(plt, 400, 400, 50, 450)
#
# fig1, ax1 = plt.subplots()
# im1 = plt.imshow(binarydata_close, cmap='Greys', interpolation='nearest')
# fp.format_plot(plt, 400, 400, 450, 450)
#
# fig1, ax1 = plt.subplots()
# im1 = plt.imshow(doubleskeleton, cmap='Greys', interpolation='nearest')
# fp.format_plot(plt, 400, 400, 850, 50)

plt.show()
