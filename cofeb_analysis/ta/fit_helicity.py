# @Author: Jenkins Alec <alec>
# @Date:   2017-07-10T17:48:10-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-29T10:41:25-07:00



import numpy as np
import matplotlib.pyplot as plt
import json

import stray_field_calc_fast as sfcf
import linecut

# def fit_helicity(scannum, cropsize):

scannum = 1760
cropsize = 11

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
scan_path = path+str(scannum)+'/'
scan_params_path = scan_path+'/'+'scan_parameters.json'
material_params_path = path+'material_parameters.json'
with open(scan_params_path, 'r') as fread:
    scan_params = json.load(fread)
with open(material_params_path, 'r') as fread:
    material_params = json.load(fread)

Ms = material_params['Ms']
t = material_params['t']
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


bNV = np.loadtxt(scan_path+'bNV.txt', delimiter=',')
bNVcrop = bNV[ycenter-cropsize:ycenter+cropsize+1, xcenter-cropsize:xcenter+cropsize+1]

mz = np.loadtxt(scan_path+'mz.txt', delimiter=',')
phi_perp = np.loadtxt(scan_path+'phi_perp.txt', delimiter=',')

gridsize = len(mz)
resRatio = int(gridsize/xres)


helicity_list = np.arange(60,120,2)*np.pi/180
# helicity_list = [70*np.pi/180]

hnum = len(helicity_list)
helicity_chisq = np.zeros((hnum))

for i in range(len(helicity_list)):
    print('getting chi^2 for helicty angle = '+str(helicity_list[i]*180/np.pi))
    mxh = -np.sqrt(1-mz**2)*np.cos(phi_perp + helicity_list[i])
    myh = -np.sqrt(1-mz**2)*np.sin(phi_perp + helicity_list[i])
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

mxh = -np.sqrt(1-mz**2)*np.cos(phi_perp + helicity_list[hminarg])
myh = -np.sqrt(1-mz**2)*np.sin(phi_perp + helicity_list[hminarg])
scd, vcd, meff, hk, h = sfcf.stray_field_calc_fast(mxh,myh,mz,Ms*t,scanSize,height)

h[2] = h[2] + zfield

savepath = scan_path+'/stray_field_sim/'
np.savetxt(savepath+'helicity_x_'+str(scannum)+'.txt', h[0], delimiter=',')
np.savetxt(savepath+'helicity_y_'+str(scannum)+'.txt', h[1], delimiter=',')
np.savetxt(savepath+'helicity_z_'+str(scannum)+'.txt', h[2], delimiter=',')

bNVsim = np.abs(h[0]*np.sin(theta)*np.cos(phi) + h[1]*np.sin(theta)*np.sin(phi)
                + (h[2]+zfield)*np.cos(theta))
bNVsim = bNVsim*(1e4) # convert to G
bNVsimLowRes = bNVsim[::resRatio, ::resRatio]

cutSize = 2.2
scanSize =scanSize*(1e6)

fieldCutNV = linecut.linecut(bNVsim, scanSize, cutSize, 0)
bNVCut = linecut.linecut(bNV, scanSize, cutSize, 0, xcenter, ycenter)

fig1, ax1 = plt.subplots()
plt.plot(fieldCutNV[0], fieldCutNV[1], label=r'min $\chi^2$ helicity angle')
plt.plot(bNVCut[0],bNVCut[1], label=u'data')
plt.legend(loc=2,borderaxespad=1,prop={'size':10})

plt.show()
