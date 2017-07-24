# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np
from scipy import ndimage
from scipy import misc
from scipy import signal
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import json
import ctypes
import time
from PIL import Image

import load_scan as lscan
import vector_reconstruction as vr
import fourier_image as fi

pi = np.pi


# def vector_m_reconstruction(scannum):

scannum = 1760

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
scan_params_path = path+str(scannum)+'/'+'scan_parameters.json'
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

data = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,yres,Dgs=2870,maxfgrad=20,fieldangle=True)
misc.imsave('/Users/alec/UCSB/scan_images/full-field/ff'+str(scannum)+'.png', data[0])
ffmask = ndimage.imread('/Users/alec/UCSB/scan_images/full-field/ff'+str(scannum)+'mask.png',flatten = True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

np.savetxt(path+str(scannum)+'/fail_index.txt', data[7])

print(data[6])

#---------------- FIT FUNCTIONS ----------------------------------
#-----------------------------------------------------------------

def fit_tanh(x, *params):
    y = np.zeros_like(x)
    c = params[0]
    a = params[1]
    x0 = params[2]
    wid = params[3]

    y = c+(a/2)*np.tanh((x-x0)/wid)
    return y


#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

datas = np.multiply(ffmask,data[0])
datas0 = np.add(datas,-np.cos(theta)*zfield)
datas0filter = signal.wiener(datas0, 3)
datas0filter = signal.medfilt(datas0filter, 5)

wimdata = fi.window_image(datas0)
wimdatafilter = fi.window_image(datas0filter, power=1/2)

recon_data = vr.vector_reconstruction(wimdata, data[1], theta, thetaError, phi, height, scanSize, kcutoff=1)

bxdata = recon_data[0]
bydata = recon_data[1]
bzdata = recon_data[2]+zfield
mdata = recon_data[3]
Vdata = recon_data[4]
bzdataError = recon_data[5]

mdataint = ndimage.interpolation.zoom(mdata, 2, order=1)

np.savetxt(path+'bNV.txt', data[0], delimiter=',')
np.savetxt(path+'bNVError.txt', data[1], delimiter=',')
np.savetxt(path+'bz.txt', bzdata, delimiter=',')
np.savetxt(path+'bzError.txt', bzdataError, delimiter=',')
np.savetxt(path+'bx.txt', bxdata, delimiter=',')
np.savetxt(path+'by.txt', bydata, delimiter=',')
np.savetxt(path+'mdata.txt', mdata, delimiter=',')
np.savetxt(path+'V.txt', Vdata, delimiter=',')

#---------------- LINECUTS ---------------------------------------
#-----------------------------------------------------------------

phinum = 30
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

r0s = np.append(r0s,r0s[0])
phis = np.linspace(0, 2*np.pi, phinum+1)

r0s_int = interpolate.CubicSpline(phis, r0s, bc_type='periodic')

phi_int_num = 200
phis_int_list = np.linspace(min(phis), max(phis), phi_int_num)
x0_int, y0_int = (r0s_int(phis_int_list)*np.cos(phis_int_list), r0s_int(phis_int_list)*np.sin(phis_int_list))
x0_int_ex = np.append(x0_int, x0_int[0])
y0_int_ex = np.append(y0_int, y0_int[0])
dw_grad = -(y0_int_ex[1:]-y0_int_ex[0:-1]), (x0_int_ex[1:]-x0_int_ex[0:-1])
theta_grad = np.arctan2(dw_grad[1], dw_grad[0])

gridnum = 200
xm = np.linspace(-scanSize/2, scanSize/2, gridnum, endpoint=False)
ym = np.linspace(-scanSize/2, scanSize/2, gridnum, endpoint=False)
xm_grid, ym_grid = np.meshgrid(xm, ym)
mindist = np.zeros((gridnum, gridnum), dtype=np.double)
theta_grad_grid = np.zeros((gridnum, gridnum), dtype=np.double)

dw_dist_lib = ctypes.cdll.LoadLibrary('./dw_dist.so')
get_dw_dist = dw_dist_lib.dw_dist

t1 = time.time()
get_dw_dist(ctypes.c_void_p(xm_grid.ctypes.data), ctypes.c_void_p(ym_grid.ctypes.data),
            ctypes.c_void_p(x0_int.ctypes.data), ctypes.c_void_p(y0_int.ctypes.data),
            ctypes.c_void_p(theta_grad.ctypes.data), ctypes.c_int(gridnum), ctypes.c_int(phi_int_num),
            ctypes.c_void_p(mindist.ctypes.data), ctypes.c_void_p(theta_grad_grid.ctypes.data))
t2 = time.time()

print(t2-t1)

mz = np.tanh(mindist/DWWidth)

# Bloch
mxBloch = np.sqrt(1-mz**2)*np.sin(theta_grad_grid)
myBloch = np.sqrt(1-mz**2)*np.cos(theta_grad_grid)

# left-handed Neel
mxLNeel = np.sqrt(1-mz**2)*np.cos(theta_grad_grid)
myLNeel = np.sqrt(1-mz**2)*np.sin(theta_grad_grid)

# right-handed Neel
mxRNeel = -np.sqrt(1-mz**2)*np.cos(theta_grad_grid)
myRNeel = -np.sqrt(1-mz**2)*np.sin(theta_grad_grid)

savepath = path+str(scannum)+'/stray_field_sim/'
np.savetxt(savepath+'mz.txt', mz, delimiter=',')
np.savetxt(savepath+'theta_grad_grid.txt', theta_grad_grid, delimiter=',')
np.savetxt(savepath+'mxBloch.txt', mxBloch, delimiter=',')
np.savetxt(savepath+'myBloch.txt', myBloch, delimiter=',')
np.savetxt(savepath+'mxLNeel.txt', mxLNeel, delimiter=',')
np.savetxt(savepath+'myLNeel.txt', myLNeel, delimiter=',')
np.savetxt(savepath+'mxRNeel.txt', mxRNeel, delimiter=',')
np.savetxt(savepath+'myRNeel.txt', myRNeel, delimiter=',')

# plt.close('all')
#
# savepath = path+str(scannum)+'/'
#
# fig, ax = plt.subplots()
# plt.imshow(bxdata, interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# plt.savefig(savepath+'Bx_recon.pdf',  bbox_inches='tight')
#
# fig, ax = plt.subplots()
# plt.imshow(bydata, interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# plt.savefig(savepath+'By_recon.pdf',  bbox_inches='tight')
#
# fig, ax = plt.subplots()
# plt.imshow(bzdata, interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# plt.savefig(savepath+'Bz_recon.pdf',  bbox_inches='tight')
#
#
# bpdata = np.sqrt(bxdata**2 + bydata**2)
#
# fig, ax = plt.subplots()
# plt.imshow(bpdata, interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# plt.savefig(savepath+'Bp_recon.pdf',  bbox_inches='tight')
#
#
# fig, ax = plt.subplots()
# plt.imshow(data[4], interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# plt.savefig(savepath+'Bp_fit.pdf',  bbox_inches='tight')
#
#
# plt.show()
