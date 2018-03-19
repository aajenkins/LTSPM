# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math as math
import matplotlib.gridspec as gridspec
import cwesr_fit_single as cwesr
import plotting.format_plots_tkagg as fp
#import matplotlib.gridspec as gridspec

filename = 'linecut-middle-20uA-0222'
scannum = 2028
path = '/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/'+filename
fitdata = []

num_avg = 50
xaxis_scale = 1e-3
yaxis_scale = 1e-3

filelist = np.arange(107,136)

filelen = len(filelist)

filenum=107
filepath = path+str(filenum).zfill(6)
data = np.loadtxt(filepath+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]

edata = np.transpose(data)

x, y = edata

cwresult = cwesr.cwesr_fit_N15(x, y, gwidth=1.5, gctr = 2902)
popt = cwresult[0]
fit = cwresult[2]
fitg = cwresult[3]

xs = [xaxis_scale*k for k in x]
ys = [yaxis_scale*k for k in y]
fits = [yaxis_scale*k for k in fit]
fitsg = [yaxis_scale*k for k in fitg]


fig, ax = plt.subplots(figsize=(6,4))

plt.plot(xs, ys)
plt.plot(xs, fits, '-r')
plt.plot(xs, fitsg, '-b')
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.005))
plt.subplots_adjust(left=0.17, bottom=0.15, right=0.98, top=1.0, wspace=0, hspace=0)

plt.show()
