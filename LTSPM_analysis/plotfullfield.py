# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:17:32 2016

@author: alec
"""
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import load_scan as ls

font = {'family': 'Arial',
        'weight': 'normal',
        'size': 15}

matplotlib.rc('font', **font)

# file constants
dres = 60
vsize = 0.8
dsize = vsize*5
filenum = 1847

ffdata = ls.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,maxfgrad=15)

plt.close('all')

plt.figure(1, [5, 5])

imgplot = plt.imshow(ffdata[0], cmap='bone', interpolation='nearest', extent=[0, dsize, 0, dsize])

plt.colorbar(fraction=0.046, pad=0.04)
plt.show()
plt.tight_layout()
pylab.savefig('/Users/alec/UCSB/scan_images/full-field/ff_'+str(filenum)+'.png',bbox_inches='tight')

fig = plt.gcf()
fig.canvas.manager.window.raise_()
