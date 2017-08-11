# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np
from scipy import ndimage
from scipy import signal
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import ctypes
from numpy.ctypeslib import ndpointer
import time
import json

import plotting.format_plots_tkagg as fp

pi = np.pi


#---------------- FIT FUNCTION -----------------------------------
#-----------------------------------------------------------------

def fit_tanh(x, *params):
    y = np.zeros_like(x)
    c = params[0]
    a = params[1]
    x0 = params[2]
    wid = params[3]

    y = c+(a/2)*np.tanh((x-x0)/wid)
    return y


def calc_DW(scannum, simGridSize=400, dwIntSize=200, helicities=[0,90,180],
subtractMsSat=False):

    helicities = np.array(helicities)*pi/180

    path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
    scan_path = path+str(scannum)+'/'
    scan_params_path = scan_path+'scan_parameters.json'
    material_params_path = path+'material_parameters.json'
    with open(scan_params_path, 'r') as fread:
        scan_params = json.load(fread)
    with open(material_params_path, 'r') as fread:
        material_params = json.load(fread)

    Ms = material_params['Ms']
    t = material_params['t']
    DWWidth = material_params['DWWidth']

    scanSize = scan_params['scanSize']
    xres = scan_params['xres']
    yres = scan_params['yres']
    xcenter = scan_params['xcenter']
    ycenter = scan_params['ycenter']

    mdata = np.loadtxt(scan_path+'mdata.txt', delimiter=',')

    #---------------- LINECUTS ---------------------------------------
    #-----------------------------------------------------------------

    phinum = 15
    lcnum = 15
    lclen = 15
    mphi = np.zeros((phinum,lcnum))
    for i in range(0,phinum):
        phi = i*2*pi/phinum
        x1, y1 = xcenter+lclen*np.cos(phi), ycenter+lclen*np.sin(phi)
        x, y = np.linspace(xcenter, x1, lcnum), np.linspace(ycenter, y1, lcnum)
        mphi[i] = ndimage.map_coordinates(np.transpose(mdata), np.vstack((x,y)), order=1)

    #---------------- LINE FITS --------------------------------------
    #-----------------------------------------------------------------

    xf = np.arange(0, lcnum)
    fits = np.zeros((phinum, lcnum))
    guesses = np.zeros((phinum, 4))
    widths = np.zeros(phinum)
    r0s = np.zeros(phinum)
    angles = np.linspace(0, 2*pi, phinum+1)

    for i in range (0,phinum):
        y = mphi[i]
        guesses[i] = [(y[-1]+y[0])/2,y[-1]-y[0],6,1]
        popt, pcov = curve_fit(fit_tanh, xf, mphi[i], p0=guesses[i])
        fits[i] = fit_tanh(xf, *popt)
        widths[i] = np.abs(popt[3])
        r0s[i] = popt[2]*scanSize/xres

    phis = np.linspace(0, 2*np.pi, phinum, endpoint=False)

    r0s = np.append(r0s[-4:],np.append(r0s,r0s[:4]))
    phis = np.append(phis[-4:]-2*pi,np.append(phis,phis[:4]+2*pi))

    r0s_int = interpolate.UnivariateSpline(phis, r0s, s=0, k=4)
    dr0s_int = r0s_int.derivative()

    phis_int_list = np.linspace(0, 2*pi, dwIntSize, endpoint=False)
    x0_int, y0_int = (r0s_int(phis_int_list)*np.cos(phis_int_list), r0s_int(phis_int_list)*np.sin(phis_int_list))
    s0_int = np.sqrt((x0_int[1:]-x0_int[:-1])**2 + (y0_int[1:]-y0_int[:-1])**2)
    s0_int = np.append(s0_int, np.sqrt((x0_int[0]-x0_int[-1])**2 + (y0_int[0]-y0_int[-1])**2) )

    xm = np.linspace(-scanSize/2, scanSize/2, simGridSize, endpoint=False)
    ym = np.linspace(-scanSize/2, scanSize/2, simGridSize, endpoint=False)
    xm_grid, ym_grid = np.meshgrid(xm, ym)
    phi_grid = np.mod(np.arctan2(ym_grid, xm_grid),2*pi)

    mindist = np.zeros((simGridSize, simGridSize), dtype=np.double)

    dw_dist_lib = ctypes.cdll.LoadLibrary('./dw_dist.so')
    get_dw_dist = dw_dist_lib.dw_dist
    get_dw_dist.argtypes = [ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), ndpointer(ctypes.c_double),
                            ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_int,
                            ndpointer(ctypes.c_double)]


    t1 = time.time()
    get_dw_dist(np.ascontiguousarray(xm_grid, np.double), np.ascontiguousarray(ym_grid, np.double),
                np.ascontiguousarray(x0_int, np.double), np.ascontiguousarray(y0_int, np.double),
                np.ascontiguousarray(s0_int, np.double), simGridSize, dwIntSize,
                np.ascontiguousarray(mindist, np.double))

    t2 = time.time()

    print(t2-t1)

    msign = np.round( (np.sqrt(xm_grid**2 + ym_grid**2) - (r0s_int(phi_grid)))
                     / np.abs(np.sqrt(xm_grid**2 + ym_grid**2) - r0s_int(phi_grid)) )
    mz = np.tanh(msign*mindist/DWWidth)

    if (subtractMsSat):
        MsSatScanNum = scan_params['MsSatScanNum']
        scan_sat_path = path+str(MsSatScanNum)+'/'
        mzSat = np.loadtxt(scan_sat_path+'mdata.txt', delimiter=',')
        mzSatInt = ndimage.interpolation.zoom(mzSat, int(simGridSize/len(mzSat)), order=3)
        print(np.mean(mzSat))
        # mz = mz * (np.divide(mzSatInt, np.mean(mzSat)))
        mz = mz * mzSatInt

    mzgrad = np.gradient(mz)

    mzgradx = mzgrad[1]
    mzgrady = mzgrad[0]

    phi_perp = pi + np.arctan2(mzgrady, mzgradx)

    np.savetxt(scan_path+'mz.txt', mz, delimiter=',')
    np.savetxt(scan_path+'phi_perp.txt', phi_perp, delimiter=',')

    numHelcities = len(helicities)
    mxs = np.empty(numHelcities, dtype=object)
    mys = np.empty_like(mxs)
    vcd = np.empty_like(mxs)
    scd = np.empty_like(mxs)

    for i in range(numHelcities):

        print('calculating M for helicity '+str(helicities[i]))

        mxs[i] = np.sqrt(1-mz**2)*np.cos(phi_perp + helicities[i])
        mys[i] = np.sqrt(1-mz**2)*np.sin(phi_perp + helicities[i])

        vcd[i] = -(np.gradient(mxs[i])[1]+np.gradient(mys[i])[0])
        scd[i] = mz
        np.savetxt(scan_path+'mx'+str(int(np.round(helicities[i]*180/pi)))+'.txt', mxs[i], delimiter=',')
        np.savetxt(scan_path+'my'+str(int(np.round(helicities[i]*180/pi)))+'.txt', mys[i], delimiter=',')


if __name__ == "__main__":
    import sys
    if (len(sys.argv) == 2):
        calc_DW(int(sys.argv[1]))
    elif (len(sys.argv) == 3):
        calc_DW(int(sys.argv[1]), simGridSize=int(sys.argv[2]))
    elif (len(sys.argv) == 4):
        calc_DW(int(sys.argv[1]), simGridSize=int(sys.argv[2]),
                                dwIntSize=int(sys.argv[3]))
    elif (len(sys.argv) == 5):
        calc_DW(int(sys.argv[1]), simGridSize=int(sys.argv[2]),
                        dwIntSize=int(sys.argv[3]), helicities=np.array(eval(sys.argv[4])))
    elif (len(sys.argv) == 6):
        calc_DW(int(sys.argv[1]), simGridSize=int(sys.argv[2]),
                        dwIntSize=int(sys.argv[3]), helicities=np.array(eval(sys.argv[4])),
                        subtractMsSat=int(sys.argv[5]))
    else:
        print('enter scan number, mask(opt), printError(opt)')
