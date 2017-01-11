# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import misc
from scipy.optimize import curve_fit
import matplotlib.pylab as pylab
import load_scan as ls
import fourier_image as fi

pi = np.pi

scannum = 1760
#scanbacknum = 1740
xres = 50
yres = 50
zfield = 9.5
scanL = 0.6*5000
height = 66


data = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,15)

#---------------- FIND NV zero-field contour ----------------------------------
#-----------------------------------------------------------------



#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

plt.figure(1,[4.5,4.5])
plt.imshow(data[0], cmap='gray', interpolation='nearest')
plt.colorbar(fraction=0.046, pad=0.04)
plt.tight_layout()
fig = plt.gcf().canvas.manager.window
geom = fig.geometry()
x,y,dx,dy = geom.getRect()
fig.setGeometry(50,50,dx, dy)
fig.raise_()