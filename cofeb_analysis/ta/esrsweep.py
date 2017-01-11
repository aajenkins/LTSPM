# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 11:16:45 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 14}
matplotlib.rc('font', **font)

savefolder = '/Users/alec/UCSB/presentations/MRS_2016/'

for i in range(19,25):
    esrnum = str(i)
    
    datapath = '/Users/alec/UCSB/scan_data/fixed_tip_esr_field_sweep/esr0000'+esrnum+'_40.txt'
    savepath = savefolder+'esrsweep'+esrnum+'.pdf'
    esrdata = scanlist = np.loadtxt(datapath, skiprows=1, delimiter='\t')[:,[0,2]]
    esrdata = np.transpose(esrdata)
    esrdata[0] = esrdata[0]/1000
    maxfl = np.max(esrdata[1])
    esrdata[1] = esrdata[1]/maxfl
    
    plt.close('all')
    
    plt.figure(1,[3,2.5])
    
    plt.plot(esrdata[0],esrdata[1])
    plt.locator_params(nbins=4)
    plt.tight_layout()
    fig = plt.gcf()
    fig.canvas.manager.window.raise_()
    
    #plt.figure(2,[6,6])    
    #plt.hist(scandata.ravel(), bins=256, range=(0, 7e4), fc='k', ec='k')
    
    pylab.savefig(savepath,bbox_inches='tight')