# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:57:37 2016

@author: alec
"""
import numpy as np
import re


def load_ff (path, xres, yres, maxfgrad=15, maxfield=150, neighbors=4):
    
    scanlist = np.loadtxt(path, delimiter=',')[:,[2,8]]
    scanlisterr = np.loadtxt(path, delimiter=',')[:,[3,9]]
    
    scan2d = np.zeros((yres, xres))
    scan2derr = np.zeros((yres, xres))
    rscan2d = np.zeros((yres, xres))
    rscan2derr = np.zeros((yres, xres))
    
    for j in range(0,yres):
        for i in range(0,xres):
            scan2d[j,i] = abs(scanlist[i+(2*j*xres),1]-scanlist[i+(2*j*xres),0])/(2*2.8)
            scan2derr[j,i] = np.sqrt((scanlisterr[i+(2*j*xres),1]**2)+(scanlisterr[i+(2*j*xres),0]**2))
            rscan2d[j,i] = abs(scanlist[i+((2*j+1)*xres),1]-scanlist[i+((2*j+1)*xres),0])/(2*2.8)
            rscan2derr[j,i] = np.sqrt((scanlisterr[i+((2*j+1)*xres),1]**2)+(scanlisterr[i+((2*j+1)*xres),0]**2))
            if scan2d[j,i] > maxfield:
                scan2d[j,i] = scan2d[j,i-1]
                scan2derr[j,i] = scan2derr[j,i-1]
            if scan2derr[j,i] > maxfield:
                scan2d[j,i] = scan2d[j-1,i]
                scan2derr[j,i] = scan2derr[j,i-1] 
            if rscan2d[j,i] > maxfield:
                rscan2d[j,i] = rscan2d[j,i-1]
                rscan2derr[j,i] = rscan2derr[j,i-1]
            if rscan2derr[j,i] > maxfield:
                rscan2d[j,i] = rscan2d[j,i-1]
                rscan2derr[j,i] = rscan2derr[j,i-1]
    
    for j in range(0,yres):
        for i in range(0,xres):
            ni = np.array([-1,1])
            nj = np.array([-1,1])
            if (i==0):
                ni = ni[[1]]
            if (i==xres-1):
                ni = ni[[0]]
            if (j==0):
                nj=nj[[1]]
            if (j==yres-1):
                nj=nj[[0]]
            
            nilen = len(ni)
            njlen = len(nj)
            nlen = njlen+nilen
            
            mean = 0
            meanerr = 0
            rmean = 0
            rmeanerr = 0
            for k in range(njlen):
                for l in range(nilen):
                    mean = mean+(scan2d[j+nj[k],i+ni[l]]/nlen)
                    meanerr = mean+(scan2derr[j+nj[k],i+ni[l]]/nlen)
                    rmean = mean+(rscan2d[j+nj[k],i+ni[l]]/nlen)
                    rmeanerr = mean+(rscan2derr[j+nj[k],i+ni[l]]/nlen)        
                
            if abs(scan2d[j,i]-mean) > maxfgrad:
                scan2d[j,i]=mean
                scan2derr[j,i]=meanerr
            if abs(rscan2d[j,i]-rmean) > maxfgrad:
                rscan2d[j,i]=rmean
                rscan2derr[j,i]=rmeanerr
    
    return scan2d, scan2derr, rscan2d, rscan2derr

    
def load_rf_track (path, xres, yres):
    scanlist = np.loadtxt(path, delimiter='\t')[:,6]
    
    scan2d = np.zeros((yres, xres))    
    
    for i in range(0,xres):
        for j in range(0,yres):
              scan2d[j,i] = abs(2870-scanlist[(i*yres)+j])/(2.8)
    return scan2d
    
def load_contour (path):
    scanlist = np.loadtxt(path, delimiter='\t')[:,6]
    res = int(np.sqrt(len(scanlist)))
    scan2d = np.zeros((res, res))
    for j in range(0,res):
        for i in range(0,res):
            scan2d[i,j] = scanlist[i+(j*res)]
    return scan2d
    
def get_scan_size (path):
    info = open(path,'r')
    lines = info.readlines()
    size = re.search('\d+\.*\d*',lines[1])
    return float(size.group(0))*1.125
    