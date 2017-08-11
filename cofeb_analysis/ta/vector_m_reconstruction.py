# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec

imports fitted CWESR data for full-field scan, calculates
field vector components, and normalized Mz component in
Bloch gauge

input args are: scan number, bool to use a mask to assign
field polarity to the NV image, bool to print errors in the B field
calculation
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np
from scipy import ndimage
from scipy import misc
from scipy import signal
from scipy import interpolate
from scipy.optimize import curve_fit
import json
from PIL import Image

import load_scan as lscan
import vector_reconstruction as vr
import fourier_image as fi
import plotting.format_plots_tkagg as fp

pi = np.pi


def vector_m_reconstruction(scannum, kcutoff=1, mask = True, printError = False):

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

    phi = scan_params['phi']
    theta = scan_params['theta']
    thetaError = scan_params['thetaError']
    height = scan_params['height']
    heightError = scan_params['heightError']
    scanSize = scan_params['scanSize']
    xres = scan_params['xres']
    yres = scan_params['yres']
    xcenter = scan_params['xcenter']
    ycenter = scan_params['ycenter']
    zfield = scan_params['zfield']*(1e4) # convert to G

    data = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,Dgs=2870,maxfgrad=20,fieldangle=True,printNVCalcError=printError)
    misc.imsave('/Users/alec/UCSB/scan_images/full-field/ff'+str(scannum)+'.png', data[0])

    np.savetxt(scan_path+'/fail_index.txt', data[7])

    #---------------- RECONSTRUCTION ---------------------------------
    #-----------------------------------------------------------------

    datas = data[0]

    if (mask):
        ffmask = ndimage.imread('/Users/alec/UCSB/scan_images/full-field/ff'+str(scannum)+'mask.png',flatten = True)
        ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)
        datas = np.multiply(ffmask,datas)

    datas0 = np.add(datas,-np.cos(theta)*zfield)

    wimdata = fi.window_image(datas0, power=1/10)

    recon_data = vr.vector_reconstruction(wimdata, data[1], theta, thetaError, phi, height, scanSize, t, kcutoff=kcutoff)

    bxdata = recon_data[0]
    bydata = recon_data[1]
    bzdata = recon_data[2]+zfield
    mdata = recon_data[3]
    Vdata = recon_data[4]
    bzdataError = recon_data[5]

    np.savetxt(scan_path+'bNV.txt', data[0], delimiter=',')
    np.savetxt(scan_path+'bNVError.txt', data[1], delimiter=',')
    np.savetxt(scan_path+'bz.txt', bzdata, delimiter=',')
    np.savetxt(scan_path+'bzError.txt', bzdataError, delimiter=',')
    np.savetxt(scan_path+'bx.txt', bxdata, delimiter=',')
    np.savetxt(scan_path+'by.txt', bydata, delimiter=',')
    np.savetxt(scan_path+'mdata.txt', mdata, delimiter=',')
    np.savetxt(scan_path+'V.txt', Vdata, delimiter=',')


if __name__ == "__main__":
    import sys
    if (len(sys.argv) == 2):
        vector_m_reconstruction(int(sys.argv[1]))
    elif (len(sys.argv) == 3):
        vector_m_reconstruction(int(sys.argv[1]), kcutoff=float(sys.argv[2]))
    elif (len(sys.argv) == 4):
        vector_m_reconstruction(int(sys.argv[1]), kcutoff=float(sys.argv[2]), mask=eval(sys.argv[3]))
    elif (len(sys.argv) == 5):
        vector_m_reconstruction(int(sys.argv[1]), kcutoff=float(sys.argv[2]), mask=eval(sys.argv[3]),
                                printError=eval(sys.argv[4]))
    else:
        print('enter scan number, mask(opt), printError(opt)')
