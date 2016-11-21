# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:19:24 2016

@author: alec
"""

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.image as mpimg
import numpy as np
import load_scan as ls

path = '/Users/alec/UCSB/scan_data/1533/001533.scan'
res = 150
scan_size = 2.5

scandata = ls.load_contour(path,res,res)


plt.close('all')

plt.figure(1,[6,6])

imgplot = plt.imshow(scandata, cmap='gray', interpolation='nearest', extent=[-scan_size/2,scan_size/2,-scan_size/2,scan_size/2])
plt.show()
plt.tight_layout()
pylab.savefig('contour_1533.pdf',bbox_inches='tight')

fig = plt.gcf()
fig.canvas.manager.window.raise_()

