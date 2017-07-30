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
from numpy.ctypeslib import ndpointer
import time
from PIL import Image

import load_scan as lscan
import vector_reconstruction as vr
import fourier_image as fi
import format_plots_tkagg as fp

pi = np.pi


# def vector_m_reconstruction(scannum, simGridSize = 200, dwIntSize = 400, printError = False):
scannum = 1760
simGridSize = 400
dwIntSize = 200
printError = False

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
ffmask = ndimage.imread('/Users/alec/UCSB/scan_images/full-field/ff'+str(scannum)+'mask.png',flatten = True)
ffmask = np.multiply(np.add(np.multiply(ffmask,1/255),-0.5),-2)

np.savetxt(scan_path+'/fail_index.txt', data[7])

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

np.savetxt(scan_path+'bNV.txt', data[0], delimiter=',')
np.savetxt(scan_path+'bNVError.txt', data[1], delimiter=',')
np.savetxt(scan_path+'bz.txt', bzdata, delimiter=',')
np.savetxt(scan_path+'bzError.txt', bzdataError, delimiter=',')
np.savetxt(scan_path+'bx.txt', bxdata, delimiter=',')
np.savetxt(scan_path+'by.txt', bydata, delimiter=',')
np.savetxt(scan_path+'mdata.txt', mdata, delimiter=',')
np.savetxt(scan_path+'V.txt', Vdata, delimiter=',')

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
phi_perp = phi_grid - r0s_int(phi_grid)

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
mzgrad = np.gradient(mz)

# kxv = np.linspace(-pi*simGridSize/(2*scanSize), pi*simGridSize/(2*scanSize), simGridSize, endpoint=False)
# kyv = np.linspace(-pi*simGridSize/(2*scanSize), pi*simGridSize/(2*scanSize), simGridSize, endpoint=False)
# kx, ky = np.meshgrid(kxv, kyv)
# 
# k = np.sqrt(kx**2 + ky**2)
#
# fmz = np.fft.fftshift(np.fft.fft2(mz, norm="ortho"))
# mzgradx = np.real(np.fft.ifft2(np.fft.ifftshift(-1j*kx*fmz), norm="ortho"))
# mzgrady = np.real(np.fft.ifft2(np.fft.ifftshift(-1j*ky*fmz), norm="ortho"))

phi_perp = pi + np.arctan2(mzgrad[0], mzgrad[1])
# phi_perp = pi + np.arctan2(mzgrady, mzgradx)

# Bloch
mxBloch = -np.sqrt(1-mz**2)*np.sin(phi_perp)
myBloch = np.sqrt(1-mz**2)*np.cos(phi_perp)

# left-handed Neel
mxLNeel = np.sqrt(1-mz**2)*np.cos(phi_perp)
myLNeel = np.sqrt(1-mz**2)*np.sin(phi_perp)

# right-handed Neel
mxRNeel = -np.sqrt(1-mz**2)*np.cos(phi_perp)
myRNeel = -np.sqrt(1-mz**2)*np.sin(phi_perp)

vcdBloch = -(np.gradient(mxBloch)[1]+np.gradient(myBloch)[0])
vcdLNeel = -(np.gradient(mxLNeel)[1]+np.gradient(myLNeel)[0])

np.savetxt(scan_path+'mz.txt', mz, delimiter=',')
np.savetxt(scan_path+'phi_perp.txt', phi_perp, delimiter=',')
np.savetxt(scan_path+'mxBloch.txt', mxBloch, delimiter=',')
np.savetxt(scan_path+'myBloch.txt', myBloch, delimiter=',')
np.savetxt(scan_path+'mxLNeel.txt', mxLNeel, delimiter=',')
np.savetxt(scan_path+'myLNeel.txt', myLNeel, delimiter=',')
np.savetxt(scan_path+'mxRNeel.txt', mxRNeel, delimiter=',')
np.savetxt(scan_path+'myRNeel.txt', myRNeel, delimiter=',')

plt.close('all')

fig, ax = plt.subplots()
plt.imshow(mz, interpolation="Nearest")
plt.colorbar()
plt.xticks([])
plt.yticks([])

fig, ax = plt.subplots()
plt.imshow(vcdBloch, interpolation="Nearest")
plt.colorbar()
plt.xticks([])
plt.yticks([])

# fig, ax = plt.subplots()
# plt.imshow(bxdata, interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# # plt.savefig(scan_path+'Bx_recon.pdf',  bbox_inches='tight')
#
# fig, ax = plt.subplots()
# plt.imshow(bydata, interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# # plt.savefig(scan_path+'By_recon.pdf',  bbox_inches='tight')
#
# fig, ax = plt.subplots()
# plt.imshow(bzdata, interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# # plt.savefig(scan_path+'Bz_recon.pdf',  bbox_inches='tight')
#
#
# bpdata = np.sqrt(bxdata**2 + bydata**2)
#
# fig, ax = plt.subplots()
# plt.imshow(bpdata, interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# # plt.savefig(scan_path+'Bp_recon.pdf',  bbox_inches='tight')
#
#
# fig, ax = plt.subplots()
# plt.imshow(data[4], interpolation="Nearest")
# plt.colorbar()
# plt.xticks([])
# plt.yticks([])
# plt.savefig(scan_path+'Bp_fit.pdf',  bbox_inches='tight')
fp.format_plots(plt, small=False)
plt.show()

# if __name__ == "__main__":
#     import sys
#     if (len(sys.argv) == 2):
#         vector_m_reconstruction(int(sys.argv[1]))
#     elif (len(sys.argv) == 3):
#         vector_m_reconstruction(int(sys.argv[1]),simGridSize=int(sys.argv[2]))
#     elif (len(sys.argv) == 4):
#         vector_m_reconstruction(int(sys.argv[1]),simGridSize=int(sys.argv[2]),dwIntSize=int(sys.argv[3]))
#     else:
#         print('enter scan number')
