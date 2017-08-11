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

import plotting.format_plots_tkagg as fp

savefolder = '/Users/alec/UCSB/papers/tacofeb/figures/'

esrnum = str(24)

datapath = '/Users/alec/UCSB/scan_data/fixed_tip_esr_field_sweep/esr0000'+esrnum+'_40.txt'
savepath = savefolder+'esrsweep'+esrnum+'.pdf'
esrdata = scanlist = np.loadtxt(datapath, skiprows=1, delimiter='\t')[10:100,[0,2]]
esrdata = np.transpose(esrdata)
esrdata[0] = esrdata[0]/1000
maxfl = np.max(esrdata[1])
esrdata[1] = esrdata[1]/maxfl

plt.close('all')

fig, ax = plt.subplots(figsize=[3.1,2.2])

plt.plot(esrdata[0],esrdata[1],linewidth=2.0)
ax.set_xticks(np.linspace(2.83,2.91,5))
ax.set_yticks(np.linspace(0.88,1.0,4))
plt.xlim([2.82,2.92])
ax.yaxis.tick_right()
plt.xlabel(r'$2\gamma B_{NV}$')
# fp.format_plot(plt, 400, 400, 0, 50, tight=False)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.subplots_adjust(left=0.05, bottom=0.22, right=0.85, top=1.0)

plt.show()

#plt.figure(2,[6,6])
#plt.hist(scandata.ravel(), bins=256, range=(0, 7e4), fc='k', ec='k')

# pylab.savefig(savepath,bbox_inches='tight')
