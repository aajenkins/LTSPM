# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:26:44 2016

@author: alec
"""

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import cwesr_fit_single as cwesr
import plotting.format_plots_tkagg as fp
import calc_NV_field as cNV

def cwesr_plot_single_fit(filePath,filenum):

    plotdata = np.loadtxt(filePath, skiprows=1)[:,0:3:2]

    eplotdata = np.transpose(plotdata)
    x, y = eplotdata

    dlen = len(y)

    b, a = signal.butter(1, 0.5, btype='lowpass')
    yfilt = signal.filtfilt(b, a, y)

    cwresult = cwesr.cwesr_fit(x,y,filenum=filenum,d_gsplit=20,gauss=True,min_width=1,
                               max_width=30, max_splitting=400, max_ctr=2890)
    popt = cwresult[0]
    perr = np.diag(cwresult[1])
    fit = cwresult[2]
    fitg = cwresult[3]
    indexes = cwresult[4]

    f1=popt[1]+(popt[2]/2)
    f2=popt[1]-(popt[2]/2)

    b = cNV.calc_NV_field_angle(f1,f2)

    ss_res = np.sum((y - fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot)

    print(popt)

    fig1, ax1 = plt.subplots()
    plt.plot(indexes[0], indexes[1],'r.')
    plt.plot(x, y)
    # plt.plot(x, yfilt)
    plt.plot(x, fitg, '-k')
    plt.plot(x, fit, '-r')
    # plt.subplots_adjust(bottom=0.5)
    fp.format_plot(plt,fheight=550,x=500,y=22)
    plt.show()
