# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""
import numpy as np
#import matplotlib.pyplot as plt
import cwesr_fit_single as cwesr
#import matplotlib.gridspec as gridspec


def cwesr_fit_array(scannum, name, filestart, fileend, num_avg, d_gsplit=20, maxCtr=2875,
                    maxWidth=15, maxSplitting=200, gaussFit=False):
    base = '/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/'
    basepath = base+name
    savepath = base+'fitdata.txt'
    fitlogpath = base+'fitlog.txt'

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

        cwresult = cwesr.cwesr_fit(x, y, d_gsplit=d_gsplit, filenum=j, max_width=maxWidth,
                                   max_splitting=maxSplitting, max_ctr=maxCtr, gauss=gaussFit)
        popt = cwresult[0]
        pcov = cwresult[1]
        perr = np.sqrt(np.abs(np.diag(pcov)))
        poptwrite = np.append(np.append(popt[0],[popt[1]+popt[2]/2, popt[1]-popt[2]/2]), popt[3:])
        fitdata.append([poptwrite, perr])

    len_data = len(fitdata)
    print(len_data)
    f = open(savepath, 'w')
    for i in range (0, len_data):
    #    len_params = len(fitdata[i])
        for k in range(6):
            f.write(str(fitdata[i][0][k])+',')
            f.write(str(fitdata[i][1][k])+',')
        f.write(str(fitdata[i][0][6])+','+str(fitdata[i][1][6])+'\n')
    f.close()

    len_log = len(fitlog)
    ff = open(fitlogpath, 'w')
    ff.write('filestart = '+str(fitlog[0])+', fileend = '+str(fitlog[1])+'\n')
    ff.write('\n fit failures by filenum \n')
    for i in range(2, len_log):
        ff.write(str(fitlog[i])+'\n')
    ff.close()
