# @Author: Jenkins Alec <alec>
# @Date:   2017-07-10T17:48:10-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-11T19:51:16-07:00



import numpy as np
import stray_field_calc_fast as sfcf
import json

# def fit_helicity(scannum, cropsize):

scannum = 1760
cropsize = 11

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
scan_params_path = path+str(scannum)+'/'+'scan_parameters.json'
material_params_path = path+'material_parameters.json'
with open(scan_params_path, 'r') as fread:
    scan_params = json.load(fread)
with open(material_params_path, 'r') as fread:
    material_params = json.load(fread)

Ms = material_params['Ms2']
t = material_params['t2']
MstError = material_params['MstError']

theta = scan_params['theta']
phi = scan_params['phi']
height = scan_params['height']
heightError = scan_params['heightError']
scanSize = scan_params['scanSize']
xres = scan_params['xres']
yres = scan_params['yres']
xcenter = scan_params['xcenter']
ycenter = scan_params['ycenter']
zfield = scan_params['zfield']


bNV = np.loadtxt(path+'bNV.txt', delimiter=',')
bNVcrop = bNV[ycenter-cropsize:ycenter+cropsize+1, xcenter-cropsize:xcenter+cropsize+1]

mz = np.loadtxt(path+'mz.txt', delimiter=',')
theta_grad_grid = np.loadtxt(path+'theta_grad_grid.txt', delimiter=',')

gridsize = len(mz)
resRatio = int(gridsize/xres)


# helicity_list = np.arange(60,85,2)*np.pi/180
helicity_list = [70*np.pi/180]

hnum = len(helicity_list)
helicity_chisq = np.zeros((hnum))

for i in range(len(helicity_list)):
    print('getting chi^2 for helicty angle = '+str(helicity_list[i]*180/np.pi))
    mxh = -np.sqrt(1-mz**2)*np.cos(theta_grad_grid + helicity_list[i])
    myh = -np.sqrt(1-mz**2)*np.sin(theta_grad_grid + helicity_list[i])
    scd, vcd, meff, hk, h = sfcf.stray_field_calc_fast(mxh,myh,mz,Ms*t,scanSize,height)

    bNVsim = np.abs(h[0]*np.sin(theta)*np.cos(phi) + h[1]*np.sin(theta)*np.sin(phi)
                    + (h[2]+zfield)*np.cos(theta))
    bNVsim = bNVsim*(1e4) # convert to G
    bNVsimLowRes = bNVsim[::resRatio, ::resRatio]
    xscenter = int(xres/2)
    bNVsimLowRes = bNVsimLowRes[xscenter-cropsize:xscenter+cropsize+1,
                                xscenter-cropsize:xscenter+cropsize+1]
    helicity_chisq[i] = np.sum((bNVsimLowRes-bNVcrop)**2)/((xres-2*cropsize+1)**2)
    print(helicity_chisq[i])

hminarg = np.argmin(helicity_chisq)

mxh = np.sqrt(1-mz**2)*np.cos(theta_grad_grid + helicity_list[hminarg])
myh = np.sqrt(1-mz**2)*np.sin(theta_grad_grid + helicity_list[hminarg])
scd, vcd, meff, hk, h = sfcf.stray_field_calc_fast(mxh,myh,mz,Ms*t,scanSize,height)

bNVsim = np.abs(h[0]*np.sin(theta)*np.cos(phi) + h[1]*np.sin(theta)*np.sin(phi)
                + (h[2]+zfield)*np.cos(theta))
bNVsim = bNVsim*(1e4) # convert to G

savepath = '/Users/alec/UCSB/cofeb_analysis_data/ta/'+str(scannum)+'/stray_field_sim/'
