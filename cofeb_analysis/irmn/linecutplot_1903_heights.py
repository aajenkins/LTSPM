# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:17:32 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MaxNLocator
from scipy import ndimage
from scipy import misc
import numpy as np
import json
import glob

import load_scan as lscan
import calc_NV_field as cNV
import get_bnv_theta_error as gbnve
import format_plot as fp

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)

pi = np.pi

#material parameters

path = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
filespec = 'Msfixed'
cal_params_path = path+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

phi = cal_params['phi']
theta = cal_params['theta']
thetaError = cal_params['thetaError']

# surfmag = 1.068e6
sfieldpre = 1e4
zfield = 0
zfields = zfield/sfieldpre

#file constants
hnum = 1
rnum = 1
dres = 50
vsize = 0.5
dsize = vsize*5
filenum = 1903

datapath = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
simpath = datapath+'stray_field_sim/heights/'
errnames = ["lower", "mean", "upper"]
filespec = "Msfixed"

heights = np.arange(135,140,2)
numheights = len(heights)
nhfms = []
nh = []
nhcut = []
for i in range(0, numheights):
    nhfms.append([])
    nh.append([])
    nhcut.append([])

for j in range(0, numheights):
    nhfms[j] = np.append(nhfms[j],glob.glob(simpath+'h_*'+str(heights[j])+'*lowres_1903.txt'))
    nh[j] = [[],[],[]]
    for k in range(0, 3):
        nh[j][k] = sfieldpre * np.loadtxt(nhfms[j][k], delimiter=',')

#simulation constants
ssize = 2.5
slen = len(nh[0][0][0])
sres = ssize/slen

nhnv = np.zeros((numheights, slen, slen))

for m in range(0, numheights):
    for j in range(0, slen):
        for i in range(0, slen):
            nhnv[m][j,i] = cNV.calc_NV_field(nh[m][0][j,i], nh[m][1][j,i], nh[m][2][j,i], theta, phi)

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,maxfgrad=30)
bxdata = np.loadtxt(datapath+'bx_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bydata = np.loadtxt(datapath+'by_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bzdata = np.loadtxt(datapath+'bz_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bzdataError = np.loadtxt(datapath+'bzError_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')

ycenter = 28
xcenter = 22

ffdatacrop = ffdata[0][ycenter-20:ycenter+21,xcenter-20:xcenter+21]
bzdatacrop = bzdata[ycenter-20:ycenter+21,xcenter-20:xcenter+21]

nhnvcrop = [nhnv[i][5:46,5:46] for i in range(0,len(nhnv))]
nhzcrop = [nh[i][2][5:46,5:46] for i in range(0,len(nh))]

cutlength = 2.0

dcutnum = dres*cutlength/dsize
scutnum = cutlength/sres

cutcrop = [xcenter-int(np.ceil(dcutnum/2)), xcenter+int(np.floor(dcutnum/2))]

ffxcut = [np.arange((cutcrop[0]-xcenter)*dsize/dres, (cutcrop[1]-xcenter)*dsize/dres, dsize/dres),ffdata[0][ycenter,cutcrop[0]:cutcrop[1]],ffdata[1][ycenter,cutcrop[0]:cutcrop[1]]]
ffycut = [np.arange((cutcrop[0]-ycenter)*dsize/dres, (cutcrop[1]-ycenter)*dsize/dres, dsize/dres),ffdata[0][cutcrop[0]:cutcrop[1],xcenter],ffdata[1][cutcrop[0]:cutcrop[1],xcenter]]

ffcut = ffxcut
x0, x1, y0, y1 = cutcrop[0], cutcrop[1], ycenter, ycenter

sycenter = int(-1+slen/2)
sxcenter = int(slen/2)

scutcrop = [sxcenter-int(np.ceil(scutnum/2)),sxcenter+int(np.floor(scutnum/2))]

xs = np.add(np.multiply(np.arange(-ssize/2,ssize/2,sres),1),0)

for m in range(0, numheights):
    nhcut[m] = [xs,nhnv[m][sycenter,:]]

sx0, sx1, sy0, sy1 = scutcrop[0], scutcrop[1], sycenter, sycenter

# ---------------------- BZ LINECUTS -------------------------------------------------
# ------------------------------------------------------------------------------------

bx0, by0 = xcenter, ycenter
simbx0, simby0 = slen/2, slen/2
phinum = 8
bzlclen = 15
simbzlclen = int(bzlclen*dsize/(dres*sres))
bzs = np.arange(0, bzlclen*dsize/dres-0.000001, dsize/dres)
simbzs = np.arange(0, simbzlclen*sres, sres)
bzphi = np.zeros((phinum, bzlclen))
bzphiError = np.zeros((phinum, bzlclen))
simbznhphi = np.zeros((phinum, simbzlclen))

for i in range(0,phinum):
   phic = i*pi/phinum
   bx1, by1 = bx0-bzlclen*np.cos(phic), by0-bzlclen*np.sin(phic)
   bx2, by2 = bx0+bzlclen*np.cos(phic), by0+bzlclen*np.sin(phic)
   bx, by = np.linspace(bx1, bx2, bzlclen), np.linspace(by1, by2, bzlclen)
   bzphi[i] = ndimage.map_coordinates(np.transpose(bzdata), np.vstack((bx,by)), order=1)
   bzphiError[i] = ndimage.map_coordinates(np.transpose(bzdataError), np.vstack((bx,by)), order=1)
   simbx1, simby1 = simbx0-simbzlclen*np.cos(phic), simby0-simbzlclen*np.sin(phic)
   simbx2, simby2 = simbx0+simbzlclen*np.cos(phic), simby0+simbzlclen*np.sin(phic)
   simbx, simby = np.linspace(simbx1, simbx2, simbzlclen), np.linspace(simby1, simby2, simbzlclen)
   simbznhphi[i] = ndimage.map_coordinates(np.transpose(nh[0][2]), np.vstack((simbx,simby)), order=1)

# ---------------------- PLOTS -------------------------------------------------------
# ------------------------------------------------------------------------------------

savepath = '/Users/alec/UCSB/papers/irmn/figures/'
my_dpi = 96
size = 800
size_inches = size/my_dpi

plt.close('all')

fig1, ax1 = plt.subplots()
fig1.set_size_inches(6, 6)
ax1.set_axis_off()
im1 = plt.imshow(nhnv[0], cmap='bone',interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
# fp.format_plot(plt, 450, 450, 0, 50)
pylab.savefig(savepath+'nv_sim_'+str(filenum)+filespec+'.png', format='png', dpi=my_dpi)


fig1, ax1 = plt.subplots()
fig1.set_size_inches(6, 6)
ax1.set_axis_off()
im1 = plt.imshow(ffdata[0], cmap='bone',interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
pylab.savefig(savepath+'nv_'+str(filenum)+filespec+'.png', format='png', dpi=my_dpi)


fig1, ax1 = plt.subplots()
fig1.set_size_inches(6, 6)
ax1.set_axis_off()
im1 = plt.imshow(bzdata, cmap='bone',interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
pylab.savefig(savepath+'bz_'+str(filenum)+filespec+'.png', format='png', dpi=my_dpi)


fig1, ax1 = plt.subplots()
fig1.set_size_inches(6, 6)
ax1.set_axis_off()
im1 = plt.imshow(nh[0][2], cmap='bone',interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
pylab.savefig(savepath+'bz_sim_'+str(filenum)+filespec+'.png', format='png', dpi=my_dpi)


fig1, ax1 = plt.subplots()
fig1.set_size_inches(6, 4)
for i in range(0,numheights):
    plt.plot(nhcut[i][0],nhcut[i][1],linewidth=2.0)
# plt.plot(nlcut[i][0],nlcut[i][1],linewidth=2.0, color='#2D7DD2', label=u'left-handed Néel')
plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="data")
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
plt.ylim([-5,40])
plt.xlim([-1.0,1.0])
# fp.format_plot(plt, 600, 450, 0, 450)
pylab.savefig(savepath+'linecut_'+str(filenum)+filespec+'.png', dpi=my_dpi)

fig, axes = plt.subplots(ncols=2, nrows=int(phinum/2), sharex=True, sharey=True)
fig.set_size_inches(8, 5)
xlocator=MaxNLocator(prune='both', nbins=6)
ylocator=MaxNLocator(prune='both', nbins=6)
for j in range(0,2):
    for i in range(0,int(phinum/2)):
        axes[i,j].plot(simbzs, simbznhphi[int(i+(phinum/2)*j)], color='#2D7DD2', linewidth=2.0, label=u'DMI dixed helicity', zorder=-5)
        (_, caps, _) = axes[i,j].errorbar(bzs, bzphi[int(i+(phinum/2)*j)], bzphiError[int(i+(phinum/2)*j)],color='#ED1035',linewidth=1.0,fmt='.',label=r"reconstructed B$_z$", zorder=-4)
        for cap in caps:
             cap.set_markeredgewidth(1)
        axes[i,j].yaxis.set_major_locator(ylocator)
        axes[i,j].xaxis.set_major_locator(xlocator)
pylab.ylim([-40,10])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':12})
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.7, top=0.95, wspace=0, hspace=0)
pylab.savefig(savepath+'bzcuts_'+str(filenum)+filespec+'.pdf', format='pdf')


plt.show()
