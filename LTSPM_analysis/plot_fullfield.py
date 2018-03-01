# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:19:24 2016

@author: alec
"""

# import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import json

import load_scan as lscan
import plotting.format_plots_tkagg as fp
import cwesr_plot_single_fit as cpsf

scannum = 1597

scan_params_path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'+str(scannum)+'/scan_parameters.json'
with open(scan_params_path, 'r') as fread:
    scan_params = json.load(fread)

fileName = scan_params['fileName']
numAvg = scan_params['numAvg']
startFileNum = scan_params['startFileNum']
res = scan_params['xres']

scanFolder = '/Users/alec/UCSB/scan_data/'
dataPath = scanFolder+str(scannum)+'-esrdata/'+fileName


def plot_point_fit(event):
    if event.dblclick:
        xp = int(np.round(event.xdata))
        yp = int(np.round(event.ydata))
        filenum = yp*res*2 + xp + startFileNum
        filePath = dataPath+str(filenum).zfill(6)+'_'+str(numAvg)+'.txt'
        cpsf.cwesr_plot_single_fit(filePath,filenum)
        print(filenum)

data = lscan.load_ff(scanFolder+str(scannum)+'-esrdata/fitdata.txt',res,res,fieldangle=True,
                     printNVCalcError=True)

plt.close('all')

fig, ax = plt.subplots()
im = plt.imshow(data[0])
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
fp.format_plot(plt,fwidth=500,fheight=550)
plt.show()

cid = fig.canvas.mpl_connect('button_press_event', plot_point_fit)
