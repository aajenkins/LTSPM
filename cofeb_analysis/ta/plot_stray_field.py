# @Author: Jenkins Alec <alec>
# @Date:   2017-07-09T16:46:45-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-29T22:50:25-07:00


import numpy as np
import matplotlib.pyplot as plt
import json

import format_plots_tkagg as fp
from PIL import Image
import calc_NV_field as cNVf
import load_scan as lscan
import linecut

pi = np.pi

# def plot_stray_field(scannum):
scannum = 1760

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
scan_params_path = path+str(scannum)+'/'+'scan_parameters.json'
field_path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'+str(scannum)+'/stray_field_sim/'
datapath = '/Users/alec/UCSB/cofeb_analysis_data/ta/1760/'

with open(scan_params_path, 'r') as fread:
    scan_params = json.load(fread)

phi = scan_params['phi']
theta = scan_params['theta']
xres = scan_params['xres']
xcenter = scan_params['xcenter']
ycenter = scan_params['ycenter']
scanSize = (1e6)*scan_params['scanSize'] # convert to microns

dwtypes = ["RNeel", "Bloch", "LNeel"]
errnames = ["lower", "mean", "upper"]

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,xres,maxfgrad=20)

scd = np.empty((len(dwtypes)), dtype=object)
vcd = np.empty_like(scd)
meff = np.empty_like(scd)
bx = np.empty((len(dwtypes), len(errnames)), dtype=object)
by = np.empty_like(bx)
bz = np.empty_like(bx)
bNV = np.empty_like(bx)
fieldCutNV = np.empty_like(bx)

for j in range(len(dwtypes)):
    scd[j] = np.loadtxt(field_path+'scd_'+dwtypes[j]+str(scannum)+'.txt', delimiter=',')
    vcd[j] = np.loadtxt(field_path+'vcd_'+dwtypes[j]+str(scannum)+'.txt', delimiter=',')
    meff[j] = np.loadtxt(field_path+'meff_'+dwtypes[j]+str(scannum)+'.txt', delimiter=',')
    for i in range(len(errnames)):
        bx[j,i] = (1e4)*np.loadtxt(field_path+dwtypes[j]+'_x_'+errnames[i]+'_'+str(scannum)+'.txt', delimiter=',')
        by[j,i] = (1e4)*np.loadtxt(field_path+dwtypes[j]+'_y_'+errnames[i]+'_'+str(scannum)+'.txt', delimiter=',')
        bz[j,i] = (1e4)*np.loadtxt(field_path+dwtypes[j]+'_z_'+errnames[i]+'_'+str(scannum)+'.txt', delimiter=',')
        bNV[j,i] = cNVf.calc_NV_field(bx[j,i], by[j,i], bz[j,i], theta, phi)

slen = len(bx[0,0])

cutSize = 2.2
phinum = 8
philist = np.linspace(0, pi, phinum, endpoint=False)

ffdataCut = np.empty((phinum), dtype=object)
ffdataCutErr = np.empty((phinum), dtype=object)
fieldCutNV = np.empty((phinum, len(dwtypes), len(errnames)), dtype=object)

for k in range(phinum):
    ffdataCut[k] = linecut.linecut(ffdata[0], scanSize, cutSize, philist[k], xcenter, ycenter)
    ffdataCutErr[k] = linecut.linecut(ffdata[1], scanSize, cutSize, philist[k], xcenter, ycenter)
    for j in range(len(dwtypes)):
        for i in range(len(errnames)):
            fieldCutNV[k,j,i] = linecut.linecut(bNV[j,i], scanSize, cutSize, philist[k])


#---------------- PLOTS ------------------------------------------
#-----------------------------------------------------------------

plt.close('all')

fig, ax = plt.subplots(figsize=(5,5))
im = plt.imshow(bNV[0,1], interpolation='nearest')
plt.colorbar(im)

fig, ax = plt.subplots(figsize=(5,5))
im = plt.imshow(ffdata[0], interpolation='nearest')
plt.colorbar(im)

fig1, ax1 = plt.subplots()
plt.fill_between(fieldCutNV[0,0,0][0], fieldCutNV[0,0,0][1],fieldCutNV[0,0,2][1], alpha=0.6)
plt.plot(fieldCutNV[0,0,1][0], fieldCutNV[0,0,1][1], label=u'right-handed Néel')
plt.fill_between(fieldCutNV[0,1,1][0], fieldCutNV[0,1,0][1],fieldCutNV[0,1,2][1], alpha=0.6)
plt.plot(fieldCutNV[0,1,1][0],fieldCutNV[0,1,1][1],linewidth=2.0, label=u'Bloch')
plt.fill_between(fieldCutNV[0,2,0][0], fieldCutNV[0,2,0][1],fieldCutNV[0,2,2][1], alpha=0.6)
plt.plot(fieldCutNV[0,2,1][0],fieldCutNV[0,2,1][1],linewidth=2.0, label=u'left-handed Néel')
plt.errorbar(ffdataCut[0][0],ffdataCut[0][1],yerr=ffdataCutErr[0][1],color='#000000',fmt='.',label="NV ESR data", linewidth=1.0)
plt.legend(loc=2,borderaxespad=1,prop={'size':10})

fig, axes = plt.subplots(ncols=2, nrows=int(phinum/2), sharex=True, sharey=True)
for j in range(0,2):
    for i in range(0,int(phinum/2)):
        axes[i,j].errorbar(ffdataCut[int(i+(phinum/2)*j)][0], ffdataCut[int(i+(phinum/2)*j)][1], color='#000000', fmt='.', linewidth=1.0)
        axes[i,j].plot(fieldCutNV[0,0,1][0], fieldCutNV[int(i+(phinum/2)*j),0,1][1], color='#97CC04', linewidth=2.0, label=u'right-handed Néel')
        axes[i,j].plot(fieldCutNV[0,1,1][0], fieldCutNV[int(i+(phinum/2)*j),1,1][1], color='#2D7DD2', linewidth=2.0, label="Bloch")
        axes[i,j].plot(fieldCutNV[0,2,1][0], fieldCutNV[int(i+(phinum/2)*j),2,1][1], color='#F97304', linewidth=2.0, label=u'left-handed Néel')
        axes[i,j].get_yaxis().set_visible(False)
        axes[i,j].get_xaxis().set_visible(False)
        axes[i,j].text(0.05,.8,r'$\phi$ = '+'{:d}{}{:d}'.format(int(i+(phinum/2)*j),r'$\pi$/',phinum),
            horizontalalignment='left', verticalalignment='center',
            transform=axes[i,j].transAxes, fontsize=10)
axes[int((phinum/2)-1),0].get_yaxis().set_visible(True)
axes[int((phinum/2)-1),0].get_xaxis().set_visible(True)
plt.ylim([0,40])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':10})
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.7, top=0.95, wspace=0, hspace=0)

fig, ax = plt.subplots(figsize=(5,5))
im = plt.imshow(vcd[1], interpolation='nearest')
plt.colorbar(im)

fig, ax = plt.subplots(figsize=(5,5))
im = plt.imshow(scd[1], interpolation='nearest')
plt.colorbar(im)

fp.format_plots(plt, small=True)

plt.show()

# if __name__ == "__main__":
#     import sys
#     if (len(sys.argv) == 2):
#         plot_stray_field(int(sys.argv[1]))
#     else:
#         print('enter scan number')
