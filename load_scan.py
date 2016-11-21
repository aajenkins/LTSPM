# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:57:37 2016

@author: alec
"""
import numpy as np

def load_ff (path, xres, yres, maxfgrad):
    
    scanlist = np.loadtxt(path, delimiter=',')[:,[2,8]]
    scanlisterr = np.loadtxt(path, delimiter=',')[:,[3,9]]
    
    scan2d = np.zeros((yres, xres))
    scan2derr = np.zeros((yres, xres))
    
    
    for j in range(0,yres):
        for i in range(0,xres):
              scan2d[j,i] = abs(scanlist[i+(2*j*xres),1]-scanlist[i+(2*j*xres),0])/(2*2.8)
              scan2derr[j,i] = np.sqrt((scanlisterr[i+(2*j*xres),1]**2)+(scanlisterr[i+(2*j*xres),0]**2))
              if scan2d[j,i] > 150:
                scan2d[j,i] = scan2d[j-1,i]
                scan2derr[j,i] = scan2derr[j-1,i]
              if scan2derr[j,i] > 150:
                scan2d[j,i] = scan2d[j-1,i]
                scan2derr[j,i] = scan2derr[j-1,i]  
    
    for j in range(1,yres-1):
        for i in range(1,xres-1):     
            mean = (scan2d[j-1,i-1]+scan2d[j-1,i]+scan2d[j-1,i+1]+scan2d[j,i-1]+scan2d[j,i+1]
            +scan2d[j+1,i-1]+scan2d[j+1,i]+scan2d[j+1,i+1])/8
            meanerr = (scan2derr[j-1,i-1]+scan2derr[j-1,i]+scan2d[j-1,i+1]+scan2derr[j,i-1]+scan2derr[j,i+1]
            +scan2derr[j+1,i-1]+scan2derr[j+1,i]+scan2derr[j+1,i+1])/8
            if abs(scan2d[j,i]-mean) > maxfgrad:
                scan2d[j,i]=mean
                scan2derr[j,i]=meanerr
    
    return scan2d, scan2derr
    
def load_ff_linecut (path, xres, yres, maxfgrad):
    
    scanlist = np.loadtxt(path, delimiter=',')[:,[2,8]]
    scanlisterr = np.loadtxt(path, delimiter=',')[:,[3,9]]
    
    scan2df = np.zeros((2*yres, xres))
    scan2dferr = np.zeros((2*yres, xres))
    
    
    for j in range(0,2*yres):
        for i in range(0,xres):
              scan2df[j,i] = abs(scanlist[i+(j*xres),1]-scanlist[i+(j*xres),0])/(2*2.8)
              scan2dferr[j,i] = np.sqrt((scanlisterr[i+(j*xres),1]**2)+(scanlisterr[i+(j*xres),0]**2))
              if scan2df[j,i] > 150:
                scan2df[j,i] = scan2df[j,i-1]
                scan2dferr[j,i] = scan2dferr[j,i-1]
              if scan2dferr[j,i] > 150:
                scan2df[j,i] = scan2df[j,i-1]
                scan2dferr[j,i] = scan2dferr[j,i-1]  
    
    for j in range(0,2*yres):
        for i in range(1,xres-1):     
            mean = (scan2df[j,i-1]+scan2df[j,i+1])/2
            meanerr = (scan2dferr[j,i-1]+scan2dferr[j,i+1])/2
            if abs(scan2df[j,i]-mean) > maxfgrad:
                scan2df[j,i]=mean
                scan2dferr[j,i]=meanerr
                
    return scan2df, scan2dferr
    
def load_contour (path, xres, yres):
    scanlist = np.loadtxt(path, delimiter='\t')[:,6]
    scan2d = np.zeros((xres, yres))
    for j in range(0,yres):
        for i in range(0,xres):
            scan2d[i,j] = scanlist[i+(j*xres)]
    return scan2d
    