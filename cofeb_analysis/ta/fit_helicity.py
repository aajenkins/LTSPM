# @Author: Jenkins Alec <alec>
# @Date:   2017-07-10T17:48:10-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-10T10:51:21-07:00



import numpy as np
import matplotlib.pyplot as plt
import json

import stray_field_calc_fast as sfcf
import linecut
import calc_DW
import stray_field
import calc_NV_field

# def fit_helicity(scannum, cropsize):

scannum = 1739
cropsize = 11

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
scan_path = path+str(scannum)+'/'
field_path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'+str(scannum)+'/stray_field_sim/'
scan_params_path = scan_path+'/'+'scan_parameters.json'
with open(scan_params_path, 'r') as fread:
    scan_params = json.load(fread)

theta = scan_params['theta']
phi = scan_params['phi']
xres = scan_params['xres']
yres = scan_params['yres']
xcenter = scan_params['xcenter']
ycenter = scan_params['ycenter']

xscenter = int(xres/2)

bNV = np.loadtxt(scan_path+'bNV.txt', delimiter=',')
bNVcrop = bNV[ycenter-cropsize:ycenter+cropsize+1, xcenter-cropsize:xcenter+cropsize+1]

helicity_list = np.arange(71,78,1)
# helicity_list = [70*np.pi/180]

hnum = len(helicity_list)
helicity_chisq = np.zeros((hnum))

bx = np.empty((hnum), dtype=object)
by = np.empty_like(bx)
bz = np.empty_like(bx)
bNVsim = np.empty_like(bx)
bNVsimLowRes = np.empty_like(bx)

calc_DW.calc_DW(scannum, helicities=helicity_list)
stray_field.stray_field(scannum, helicities=helicity_list, errnames=["mean"])

for i in range(hnum):
    print('getting chi^2 for helicty angle = '+str(helicity_list[i]))

    bx[i] = (1e4)*np.loadtxt(field_path+'h'+str(helicity_list[i])+'_x_mean_'+str(scannum)+'.txt', delimiter=',')
    by[i] = (1e4)*np.loadtxt(field_path+'h'+str(helicity_list[i])+'_y_mean_'+str(scannum)+'.txt', delimiter=',')
    bz[i] = (1e4)*np.loadtxt(field_path+'h'+str(helicity_list[i])+'_z_mean_'+str(scannum)+'.txt', delimiter=',')
    bNVsim[i] = calc_NV_field.calc_NV_field(bx[i], by[i], bz[i], theta, phi)

    simres = len(bNVsim[i])
    resRatio = int(simres/xres)
    bNVsimLowRes[i] = bNVsim[i][::resRatio, ::resRatio]
    bNVsimLowRes[i] = bNVsimLowRes[i][xscenter-cropsize:xscenter+cropsize+1,
                                xscenter-cropsize:xscenter+cropsize+1]
    helicity_chisq[i] = np.sum((bNVsimLowRes[i]-bNVcrop)**2)/((xres-2*cropsize+1)**2)
    print(helicity_chisq[i])
