# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:57:37 2016

@author: alec
"""
import numpy as np
import re
import calc_NV_field as cNV


def load_ff (path, xres, yres, Dgs=2870, maxfgrad=25, maxfield=150, neighbors=4):

    scanlist = np.loadtxt(path, delimiter=',')[:,[2,4]]
    scanlisterr = np.loadtxt(path, delimiter=',')[:,[3,5]]

    Bnv = np.zeros((yres, xres))
    B0 = np.zeros((yres, xres))
    theta = np.zeros((yres, xres))
    BnvErr = np.zeros((yres, xres))
    rBnv = np.zeros((yres, xres))
    rB0 = np.zeros((yres, xres))
    rtheta = np.zeros((yres, xres))
    rBnvErr = np.zeros((yres, xres))

    num_fail = 0
    fail_indices = []

    for j in range(0,yres):
        for i in range(0,xres):
            # Bnv[j,i] = np.abs(scanlist[i+(2*j*xres),1]-scanlist[i+(2*j*xres),0])/(2*2.8)
            f_ctr = scanlist[i+(2*j*xres),0]
            f_split = scanlist[i+(2*j*xres),1]
            f_splitErr = scanlisterr[i+(2*j*xres),1]
            B_list = cNV.calc_NV_field_angle(f_ctr-(f_split/2), f_ctr+(f_split/2), Dgs)
            Bnv[j,i] = B_list[0]/2.8
            B0[j,i] = B_list[1]/2.8
            theta[j,i] = B_list[2]
            BnvErr[j,i] = f_split/2.8

            # rBnv[j,i] = np.abs(scanlist[i+((2*j+1)*xres),1]-scanlist[i+((2*j+1)*xres),0])/(2*2.8)
            rf_ctr = scanlist[i+((2*j+1)*xres),0]
            rf_split = scanlist[i+((2*j+1)*xres),1]
            rf_splitErr = scanlisterr[i+((2*j+1)*xres),1]
            rB_list = cNV.calc_NV_field_angle(rf_ctr-(rf_split/2), rf_ctr+(rf_split/2), Dgs)
            rBnv[j,i] = rB_list[0]/2.8
            rB0[j,i] = rB_list[1]/2.8
            rtheta[j,i] = rB_list[2]
            rBnvErr[j,i] = rf_split/2.8

            # if(np.abs(2870-(scanlist[i+(2*j*xres),1]+scanlist[i+(2*j*xres),0])/2) > 2.5):
            if(Bnv[j,i]>500/2.8):
                Bnv[j,i] = 1000
                num_fail = num_fail+1
                fail_indices.append([j,i])
                # print(i,j,scanlist[i+(2*j*xres),1],scanlist[i+(2*j*xres),0])

    Bnv, BnvErr, rBnv, rBnvErr = interpolate_fit_fails(Bnv, BnvErr, rBnv, rBnvErr, xres, yres, maxfield, maxfgrad)

    return Bnv, BnvErr, rBnv, rBnvErr, B0, theta, num_fail, fail_indices


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

def interpolate_fit_fails(Bnv, BnvErr, rBnv, rBnvErr, xres, yres, maxfield, maxfgrad):
    for j in range(0,yres):
        for i in range(0,xres):
            if Bnv[j,i] > maxfield:
                Bnv[j,i] = Bnv[j,i-1]
                BnvErr[j,i] = BnvErr[j,i-1]
            if BnvErr[j,i] > maxfield:
                Bnv[j,i] = Bnv[j-1,i]
                BnvErr[j,i] = BnvErr[j,i-1]
            if rBnv[j,i] > maxfield:
                rBnv[j,i] = rBnv[j,i-1]
                rBnvErr[j,i] = rBnvErr[j,i-1]
            if rBnvErr[j,i] > maxfield:
                rBnv[j,i] = rBnv[j,i-1]
                rBnvErr[j,i] = rBnvErr[j,i-1]

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
                    mean = mean+(Bnv[j+nj[k],i+ni[l]]/nlen)
                    meanerr = mean+(BnvErr[j+nj[k],i+ni[l]]/nlen)
                    rmean = mean+(rBnv[j+nj[k],i+ni[l]]/nlen)
                    rmeanerr = mean+(rBnvErr[j+nj[k],i+ni[l]]/nlen)

            if abs(Bnv[j,i]-mean) > maxfgrad:
                Bnv[j,i]=mean
                BnvErr[j,i]=meanerr
            if abs(rBnv[j,i]-rmean) > maxfgrad:
                rBnv[j,i]=rmean
                rBnvErr[j,i]=rmeanerr

    return Bnv, BnvErr, rBnv, rBnvErr
