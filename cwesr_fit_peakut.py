# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""

import numpy as np
#import matplotlib.pyplot as plt
import cwesr_fit_single as cwesr
#import matplotlib.gridspec as gridspec

basepath = '/Users/alec/UCSB/scan_data/1637-esrdata/esr'
savepath = '/Users/alec/UCSB/scan_data/1637-esrdata/fitdata.txt'
fitlogpath = '/Users/alec/UCSB/scan_data/1637-esrdata/fitlog.txt'
filestart = 4025
fileend = 4774
num_avg = 10


fitdata = []
fitlog = [filestart, fileend]
#splittings_edge = []

        
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
    pcov = cwresult[1]
    perr = np.sqrt(np.diag(pcov))
    fitdata.append([popt, perr])

len_data = len(fitdata)
print(len_data)
f = open(savepath, 'w')
for i in range (0,len_data):
#    len_params = len(fitdata[i])
    for k in range (0,6):
        f.write(str(fitdata[i][0][k])+',')
        f.write(str(fitdata[i][1][k])+',')
    f.write(str(fitdata[i][0][6])+','+str(fitdata[i][1][6])+'\n')  
f.close()


len_log = len(fitlog)
ff = open(fitlogpath, 'w')
ff.write('filestart = '+str(fitlog[0])+', fileend = '+str(fitlog[1])+'\n')
ff.write('\n fit failures by filenum \n')
for i in range (2,len_log):
    ff.write(str(fitlog[i])+'\n')
ff.close()

