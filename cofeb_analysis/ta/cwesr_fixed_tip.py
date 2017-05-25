# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MaxNLocator
import format_plot as fp

#import math as math
#import matplotlib.gridspec as gridspec
#import matplotlib.gridspec as gridspec
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 13}
matplotlib.rc('font', **font)

filename = 'esr'
num_avg = 40
path = '/Users/alec/UCSB/scan_data/fixed_tip_esr_field_sweep/'+filename
savepath = '/Users/alec/UCSB/papers/tacofeb/figures/fixed_tip/'
filetype = 'pdf'

fitdata = []

filenum = 24

filepath = path+str(filenum).zfill(6)
data = np.loadtxt(filepath+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]

edata = np.transpose(data)
x, y = edata
x = x/(1.0e3)
y = y/(1.0e3)

plt.close('all')
fig, ax = plt.subplots(figsize=(3,1.5))
plt.plot(x, y, color='#2D7DD2', linewidth=2.0)
locator=MaxNLocator(prune='both', nbins=5)
ax.yaxis.set_major_locator(locator)
# ax.set_xticklabels([])
# fp.format_plot(plt, 400, 250, 0, 50)
plt.show()
plt.savefig(savepath+str(filenum)+'esr.'+filetype, format=filetype)
