# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""

import numpy as np
#import matplotlib.pyplot as plt
import cwesr_fit_single as cwesr
#import matplotlib.gridspec as gridspec

basepath = '/Users/alec/UCSB/scan_data/1567-esrdata/fittrack2'
savepath = '/Users/alec/UCSB/scan_data/1567-esrdata/fitdata.txt'
fitlogpath = '/Users/alec/UCSB/scan_data/1567-esrdata/fitlog.txt'
filestart = 5025
fileend = 7024
num_avg = 4
dwelltime = 4e-2
gamp = 7e3
gwidth = 10
lbounds1 = [0,2600,1e3,4]
ubounds1 = [1e5,3100,2e5,40]
lbounds2 = [0,2600,1e3,4,2600,1e3,4]
ubounds2 = [1e5,3100,2e5,40,3100,2e5,40]
lbounds3 = [0,2600,1e3,4,2600,1e3,4,2600,1e3,4]
ubounds3 = [1e5,3100,2e5,40,2950,2e5,40,3100,2e5,40]
lbounds4 = [0,2600,1e3,4,2600,1e3,4,2600,1e3,4,2600,1e3,4]
ubounds4 = [1e5,3100,2e5,40,3100,2e5,40,3100,2e5,40,3100,2e5,40]
maxcenshift = 10
defaultf1 = 2772
defaultf2 = 2961

fitdata = []
fitlog = [filestart, fileend]
#splittings_edge = []

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

        
for j in range (filestart,fileend+1):
    if (j%100==0):
        print(j)
    filenum=j
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
    fitdata.append(popt)

len_data = len(fitdata)
print(len_data)
f = open(savepath, 'w')
for i in range (0,len_data):
#    len_params = len(fitdata[i])
    for k in range (0,6):
        f.write(str(fitdata[i][k])+',')
    f.write(str(fitdata[i][6])+'\n')  
f.close()


len_log = len(fitlog)
ff = open(fitlogpath, 'w')
ff.write('filestart = '+str(fitlog[0])+', fileend = '+str(fitlog[1])+'\n')
ff.write('\n fit failures by filenum \n')
for i in range (2,len_log):
    ff.write(str(fitlog[i])+'\n')
ff.close()

