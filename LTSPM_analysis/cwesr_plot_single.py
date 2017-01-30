# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""

import numpy as np
import matplotlib.pyplot as plt
import cwesr_fit_single as cwesr
from pylab import get_current_fig_manager
import format_plot as fp

#import math as math
#import matplotlib.gridspec as gridspec
#import matplotlib.gridspec as gridspec

scannum = 1903
filename = 'ff4'
num_avg = 2
path = '/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/'+filename

fitdata = []

filenum = 28*50*2+10

filepath = path+str(filenum).zfill(6)
data = np.loadtxt(filepath+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]

edata = np.transpose(data)
x, y = edata

cwresult = cwesr.cwesr_fit(x,y)
popt = cwresult[0]
perr = np.diag(cwresult[1])
fit = cwresult[2]
fitg = cwresult[3]
indexes = cwresult[4]

print(popt)


plt.close(3)
plt.figure(3,[6,5])
#print(cwresult[0][4]-cwresult[0][1])
#print(perr)
plt.plot(indexes[0], indexes[1],'r.')
plt.plot(x, y)
plt.plot(x, fit, '-r')
fp.format_plot(plt, 650, 200)
