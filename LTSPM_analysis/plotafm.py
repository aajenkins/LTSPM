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

scan_num = 1927
size = 5e-5
filetype='png'
scanfolder = '/Users/alec/UCSB/scan_data/'

scan_num = str(scan_num)
datapath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.scan'
infopath = scanfolder+scan_num+'/'+scan_num.zfill(6)+'.info'

afmdata = lscan.load_contour(datapath)[0]
scan_size = lscan.get_scan_size(infopath)

slen = len(afmdata)
x = np.arange(slen)

xx, yy = np.meshgrid(x,x)

X = xx.flatten()
Y = yy.flatten()

order = 1    # 1: linear, 2: quadratic
if order == 1:
    # best-fit linear plane
    A = np.c_[X, Y, np.ones(slen**2)]
    C,_,_,_ = scipy.linalg.lstsq(A, afmdata.flatten())    # coefficients

    # evaluate it on grid
    Z = C[0]*xx + C[1]*yy + C[2]

# elif order == 2:
#     # best-fit quadratic curve
#     f = np.polynomial.polynomial.polyval2d(xx, yy, afmdata)
#
#     vander = np.polynomial.polynomial.polyvander2d(xx, yy, [1,1])
#     vander = vander.reshape((-1,vander.shape[-1]))
#     f = f.reshape((vander.shape[0],))
#     C = np.linalg.lstsq(vander, f)[0]
#
#     # evaluate it on a grid
#     Z = C[4]*xx**2. + C[5]*yy**2. + C[3]*xx*yy + C[1]*xx + C[2]*yy + C[0]


afmdataFlat = -(afmdata-Z)*4300/4

plt.close('all')

fig, ax = plt.subplots()
im = plt.imshow(afmdata, cmap='bone', interpolation='nearest')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)

fig, ax = plt.subplots()
im = plt.imshow(afmdataFlat, cmap='bone', interpolation='nearest')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)
cbar.set_label('tip height (nm)')

plt.show()
