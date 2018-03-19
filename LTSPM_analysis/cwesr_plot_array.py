# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""

import numpy as np
import matplotlib.pyplot as plt
import math as math
import matplotlib.gridspec as gridspec
import cwesr_fit_single as cwesr
import plotting.format_plots_tkagg as fp
#import matplotlib.gridspec as gridspec

filename = 'linecut-middle-20uA-0222000'
scannum = 2028
path = '/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/'+filename
filenum_path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'+str(scannum)+'/fail_index.txt'
#savepath = '/Users/alec/UCSB/scan_data/830-esrdata/fitdata.txt'
fitdata = []
#back = np.loadtxt('/Users/alec/UCSB/python_analysis/copt/fullfield_images/back.txt')

num_avg = 50
xaxis_scale = 1e-3
yaxis_scale = 1e-3

filearray = np.loadtxt(filenum_path)

filearrayX, filearrayY = np.mgrid[28:33, 26:30]

print(filearrayX)
print(filearrayY)

filelist = np.array([], dtype=int)
for i in range(len(filearrayX[0,:])):
    for j in range(len(filearrayX[:,0])):
        filelist = np.append(filelist, int(filearrayY[j][i]*50*2+(filearrayX[j][i]+1)))

filelen = len(filelist)
print(filelist)

plt.close('all')
plt.figure(1,[17,10])
gs = gridspec.GridSpec(6, 6)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)

for j in range(np.min([filelen,36])):

    filenum=filelist[j]
    filepath = path+str(filenum).zfill(6)
    data = np.loadtxt(filepath+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]

    edata = np.transpose(data)
    #print(edata)
    x, y = edata
#    y = y - back
#    bpopt = [6619.24955982, 2771.38983444, 6717.65417238, 276.51001779, 3079.50430922, 5146.73378303, 352.99559504]
#    backfit = func(x, *bpopt)
#    y = y-backfit
    #print(edata[1,2])
    cwresult = cwesr.cwesr_fit(x,y)
    popt = cwresult[0]
    fit = cwresult[2]
    fitg = cwresult[3]
    dips = cwresult[4]

    csubplot = plt.subplot(gs[(j%6),math.floor(j/6)])

#    if len(sm) >= 3:
#        fcs = popt[[1,4,7]]
#        if len(sm) >= 4:
#            fcs = popt[[1,4,7,10]]
#        fcss = sorted(fcs)
#        fcmin = fcss[len(fcss)-1]
#        fcmax = fcss[0]
#        plt.axvline(x=xaxis_scale*fcmin,color='k',ls='dashed')
#        plt.axvline(x=xaxis_scale*fcmax,color='k',ls='dashed')
    if (len(dips)>=1):
        ppxs = [xaxis_scale*k for k in dips[0]]
        ppys = [yaxis_scale*h for h in dips[1]]
        plt.plot(ppxs, ppys,'g.')
    xs = [xaxis_scale*k for k in x]
    ys = [yaxis_scale*k for k in y]
    fits = [yaxis_scale*k for k in fit]
    fitsg = [yaxis_scale*k for k in fitg]
    plt.plot(xs, ys)
    plt.plot(xs, fits, '-r')
#    plt.plot(xs, fitsg, '-b')
    start, end = csubplot.get_ylim()
    csubplot.yaxis.set_ticks(np.arange(start, end, 1))


#plt.tight_layout()
plt.show()
fig = plt.gcf()
#fig.canvas.manager.window.raise_()
#fig=plt.gcf()
#fig.canvas.manager.window.activateWindow()
#fig.canvas.manager.window.raise_()
