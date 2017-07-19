# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import cwesr_fit_single as cwesr
import format_plot as fp

#import math as math
#import matplotlib.gridspec as gridspec
#import matplotlib.gridspec as gridspec

scannum = 1760
filename = 'ff1'
num_avg = 4
path = '/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/'+filename

fitdata = []

filenum = 33#0*50*2+(0+1)

filepath = path+str(filenum).zfill(6)
plotdata = np.loadtxt(filepath+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]

eplotdata = np.transpose(plotdata)
x, y = eplotdata

dlen = len(y)

# yfilt = signal.wiener(np.append(np.append(np.flipud(y),y),np.flipud(y)), 11)
# yfilt = yfilt[dlen:2*dlen]
b, a = signal.butter(1, 0.5, btype='lowpass')
yfilt = signal.filtfilt(b, a, y)

cwresult = cwesr.cwesr_fit(x,y,filenum=filenum,d_gsplit=20)
popt = cwresult[0]
perr = np.diag(cwresult[1])
fit = cwresult[2]
fitg = cwresult[3]
indexes = cwresult[4]

print(popt)

fig1, ax1 = plt.subplots()
plt.plot(indexes[0], indexes[1],'r.')
plt.plot(x, y)
plt.plot(x, yfilt)
plt.plot(x, fitg, '-k')
plt.plot(x, fit, '-r')
plt.show()

#print(cwresult[0][4]-cwresult[0][1])
#print(perr)

# fp.format_plot(plt, 650, 200)
