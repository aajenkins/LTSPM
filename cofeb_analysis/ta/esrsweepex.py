# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 11:16:45 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
import numpy as np

import format_plot as fp

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

savefolder = '/Users/alec/UCSB/presentations/MM_2017/'

esrnum = str(24)

datapath = '/Users/alec/UCSB/scan_data/fixed_tip_esr_field_sweep/esr0000'+esrnum+'_40.txt'
savepath = savefolder+'esrsweep'+esrnum+'.pdf'
esrdata = scanlist = np.loadtxt(datapath, skiprows=1, delimiter='\t')[10:100,[0,2]]
esrdata = np.transpose(esrdata)
esrdata[0] = esrdata[0]/1000
maxfl = np.max(esrdata[1])
esrdata[1] = esrdata[1]/maxfl

plt.close('all')

plt.figure(1,[4.0,3.5])
ax1 = plt.axes()

plt.plot(esrdata[0],esrdata[1],linewidth=2.0,color='#2D7DD2')
ax1.set_xticks(np.linspace(2.82,2.92,4))
ax1.set_yticks(np.linspace(0.88,1.0,4))
plt.xlim([2.82,2.92])
# ax1.locator_params(axis='y', nbins=5)
# fp.format_plot(plt, 400, 400, 0, 50, tight=False)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

plt.show()

#plt.figure(2,[6,6])
#plt.hist(scandata.ravel(), bins=256, range=(0, 7e4), fc='k', ec='k')

pylab.savefig(savepath,bbox_inches='tight')
