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
import plotting.format_plots_tkagg as fp

#import math as math
#import matplotlib.gridspec as gridspec
#import matplotlib.gridspec as gridspec
font = {'family' : 'Verdana',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

filename = 'esr'
num_avg = 40
path = '/Users/alec/UCSB/scan_data/fixed_tip_esr_field_sweep/'+filename
savepath = '/Users/alec/UCSB/papers/tacofeb/figures/fixed_tip/'
filetype = 'pdf'

fitdata = []

filenums = [19,20,22,24]

plt.close('all')

for filenum in filenums:
    filepath = path+str(filenum).zfill(6)
    data = np.loadtxt(filepath+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]

    edata = np.transpose(data)
    x, y = edata
    x = x/(1.0e3)
    y = y/(1.0e3)

    # plt.close('all')
    fig, ax = plt.subplots(figsize=[3,1.5])
    plt.plot(x, y, linewidth=2.0)
    locator=MaxNLocator(prune='both', nbins=5, integer=True)
    ax.yaxis.set_major_locator(locator)
    ax.yaxis.tick_right()
    if (filenum != filenums[-1]):
        ax.get_xaxis().set_visible(False)
    plt.subplots_adjust(left=0.0, bottom=0.22, right=0.85, top=1.0)
    plt.show()
    plt.savefig(savepath+str(filenum)+'esr.'+filetype, format=filetype)


    # ax.set_xticklabels([])
    # fp.format_plot(plt, 400, 250, 0, 50)
