# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:19:24 2016

@author: alec
"""

# import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import scipy.linalg

import load_scan as lscan
import plotting.format_plots_tkagg as fp
from PIL import Image

scan_num = 2028
xlen = 30
size = 1.0
scanfolder = '/Users/alec/UCSB/scan_data/'

scan_num = str(scan_num)
datapath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.scan'
infopath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.info'

afmdata = np.loadtxt(datapath, delimiter='\t')[:,2]

slen = len(afmdata)
ylen = int(slen/(xlen))
x = np.arange(slen)
afmdataForward = []

for i in range(ylen):
    afmdataForward.append(afmdata[i:slen:ylen])

# print(afmdataForward)

plt.close('all')

fig, ax = plt.subplots()
im = plt.imshow(afmdataForward, cmap='bone', interpolation='nearest')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)

plt.show()
