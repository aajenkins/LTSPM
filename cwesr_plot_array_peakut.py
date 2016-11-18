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
#import matplotlib.gridspec as gridspec

basepath = '/Users/alec/UCSB/scan_data/1566-esrdata/fittrack2'
#savepath = '/Users/alec/UCSB/scan_data/830-esrdata/fitdata.txt'
fitdata = []

def func(x, *params):
        y = np.zeros_like(x)
        c = params[0]
        for i in range(1, len(params), 3):
            ctr = params[i]
            amp = params[i+1]
            wid = params[i+2]
            #y = y + amp * (np.exp( -((x - ctr-2.3)/wid)**2)+np.exp( -((x - ctr)/wid)**2)+np.exp( -((x - ctr+2.3)/wid)**2)) + c
            y = y - abs((amp * (wid/2)**2)/((x-ctr)**2+(wid/2)**2))
        y=y+c
        return y

num_avg = 3      
filestart = 2634
xaxis_scale = 1e-3
yaxis_scale = 1e-3
gamp = 10e3
gwidth = 10
lbounds1 = [0,2600,1e3,4]
ubounds1 = [1e5,3100,2e5,40]
lbounds2 = [0,2600,1e3,4,2600,1e3,4]
ubounds2 = [1e5,3100,2e5,40,3100,2e5,40]
lbounds3 = [0,2600,1e3,4,2600,1e3,4,2600,1e3,4]
ubounds3 = [1e5,3100,2e5,40,2950,2e5,40,3100,2e5,40]
lbounds4 = [0,2600,1e3,4,2600,1e3,4,2600,1e3,4,2600,1e3,4]
ubounds4 = [1e5,3100,2e5,40,3100,2e5,40,3100,2e5,40,3100,2e5,40]
maxcenshift = 40
defaultf1 = 2772
defaultf2 = 2961

plt.close('all')
plt.figure(1,[17,10])
gs = gridspec.GridSpec(5, 6)
gs.update(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.25, hspace=0.25)
        
for j in range (0,30):
   
    filenum=j+filestart
    filepath = basepath+str(filenum).zfill(6)
    data = np.loadtxt(filepath+'_'+str(num_avg)+'.txt', skiprows=1)[:,0:3:2]

    edata = np.transpose(data)
    #print(edata)
    x, y = edata
#    bpopt = [6619.24955982, 2771.38983444, 6717.65417238, 276.51001779, 3079.50430922, 5146.73378303, 352.99559504]
#    backfit = func(x, *bpopt)
#    y = y-backfit
    #print(edata[1,2])
    cwresult = cwesr.cwesr_fit(x,y)
    popt = cwresult[0]
    fit = cwresult[1]
    fitg = cwresult[2]
    dips = cwresult[3]
    
    csubplot = plt.subplot(gs[(j%5),math.floor(j/5)])

#    if len(sm) >= 3:
#        fcs = popt[[1,4,7]]
#        if len(sm) >= 4:
#            fcs = popt[[1,4,7,10]]
#        fcss = sorted(fcs)
#        fcmin = fcss[len(fcss)-1]
#        fcmax = fcss[0]
#        plt.axvline(x=xaxis_scale*fcmin,color='k',ls='dashed')  
#        plt.axvline(x=xaxis_scale*fcmax,color='k',ls='dashed')
    ppxs = [xaxis_scale*k for k in dips[0]]
    ppys = [yaxis_scale*h for h in dips[1]]
    plt.plot(ppxs, ppys,'r.')
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
#fig=plt.gcf()
#fig.canvas.manager.window.activateWindow()
#fig.canvas.manager.window.raise_()