# @Author: Jenkins Alec <alec>
# @Date:   2017-07-09T16:46:45-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-09T17:48:10-07:00


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import FormatStrFormatter

import format_plot as fp
from PIL import Image

scannum = 1760
xres = 50
yres = 50
zfield = 9.5
scanL = 0.6*(5e-6)


path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'+str(scannum)+'/'

mx = np.loadtxt(path+'mx.txt', delimiter=',')
my = np.loadtxt(path+'my.txt', delimiter=',')
mz = np.loadtxt(path+'mz.txt', delimiter=',')

#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

fig2, ax2 = plt.subplots(figsize=(5,5))
im2 = plt.imshow(mx, cmap='jet', interpolation='nearest')
plt.colorbar(im2)
fp.format_plot(plt, 400, 400, 50, 50)

fig2, ax2 = plt.subplots(figsize=(5,5))
im2 = plt.imshow(my, cmap='jet', interpolation='nearest')
plt.colorbar(im2)
fp.format_plot(plt, 400, 400, 450, 50)

fig2, ax2 = plt.subplots(figsize=(5,5))
im2 = plt.imshow(mz, cmap='jet', interpolation='nearest')
plt.colorbar(im2)
fp.format_plot(plt, 400, 400, 50, 450)

plt.show()
