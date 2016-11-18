# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:57:37 2016

@author: alec
"""
import numpy as np

def import_ff (path, xres, yres, maxfgrad):
    
    scanlist = np.loadtxt(path, delimiter=',')[:,[1,4]]
    
    scan2d = np.zeros((xres, yres))
    
    
    for j in range(0,yres):
        for i in range(0,xres):
              scan2d[j,i] = abs(scanlist[i+(2*j*50),1]-scanlist[i+(2*j*50),0])/(2*2.8)
              if scan2d[j,i] > 150:
                scan2d[j,i] = scan2d[j-1,i]
    
    for j in range(1,49):
        for i in range(1,49):     
            mean = (scan2d[j-1,i-1]+scan2d[j-1,i]+scan2d[j-1,i+1]+scan2d[j,i-1]+scan2d[j,i+1]
            +scan2d[j+1,i-1]+scan2d[j+1,i]+scan2d[j+1,i+1])/8
            if abs(scan2d[j,i]-mean) > maxfgrad:
                scan2d[j,i]=mean
    
    return scan2d