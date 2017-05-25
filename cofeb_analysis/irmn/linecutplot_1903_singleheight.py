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

font = {'family' : 'Verdana',
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

height = 140
nrightfms = glob.glob(simpath+'nr_*'+str(height)+'_1903.txt')
nleftfms = glob.glob(simpath+'nl_*'+str(height)+'_1903.txt')
blochfms = glob.glob(simpath+'b_*'+str(height)+'_1903.txt')
nlbfms = glob.glob(simpath+'h_*'+str(height)+'_1903.txt')
nrbfms = glob.glob(simpath+'hr_*'+str(height)+'_1903.txt')
nright = [[],[],[]]
nleft = [[],[],[]]
bloch = [[],[],[]]
nlb = [[],[],[]]
nrb = [[],[],[]]
for k in range(0, 3):
    nright[k] = sfieldpre * np.loadtxt(nrightfms[k], delimiter=',')
    nleft[k] = sfieldpre * np.loadtxt(nleftfms[k], delimiter=',')
    bloch[k] = sfieldpre * np.loadtxt(blochfms[k], delimiter=',')
    nlb[k] = sfieldpre * np.loadtxt(nlbfms[k], delimiter=',')
    nrb[k] = sfieldpre * np.loadtxt(nrbfms[k], delimiter=',')

#simulation constants
ssize = 2.5
slen = len(nleft[0][0])
sres = ssize/slen

nrightnv = np.zeros((slen, slen))
nleftnv = np.zeros((slen, slen))
blochnv = np.zeros((slen, slen))
nlbnv = np.zeros((slen, slen))
nrbnv = np.zeros((slen, slen))

for j in range(0, slen):
    for i in range(0, slen):
        nrightnv[j,i] = cNV.calc_NV_field(nright[0][j,i], nright[1][j,i], nright[2][j,i], theta, phi)
        nleftnv[j,i] = cNV.calc_NV_field(nleft[0][j,i], nleft[1][j,i], nleft[2][j,i], theta, phi)
        blochnv[j,i] = cNV.calc_NV_field(bloch[0][j,i], bloch[1][j,i], bloch[2][j,i], theta, phi)
        nlbnv[j,i] = cNV.calc_NV_field(nlb[0][j,i], nlb[1][j,i], nlb[2][j,i], theta, phi)
        nrbnv[j,i] = cNV.calc_NV_field(nrb[0][j,i], nrb[1][j,i], nrb[2][j,i], theta, phi)


ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,maxfgrad=30)
bxdata = np.loadtxt(datapath+'bx_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bydata = np.loadtxt(datapath+'by_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bzdata = np.loadtxt(datapath+'bz_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bzdataError = np.loadtxt(datapath+'bzError_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')

ycenter = 28
xcenter = 22

cutlength = 2.0

dcutnum = dres*cutlength/dsize
scutnum = cutlength/sres

cutcrop = [xcenter-int(np.ceil(dcutnum/2)), xcenter+int(np.floor(dcutnum/2))]

ffxcut = [np.arange((cutcrop[0]-xcenter)*dsize/dres, (cutcrop[1]-xcenter)*dsize/dres, dsize/dres),ffdata[0][ycenter,cutcrop[0]:cutcrop[1]],ffdata[1][ycenter,cutcrop[0]:cutcrop[1]]]
ffycut = [np.arange((cutcrop[0]-ycenter)*dsize/dres, (cutcrop[1]-ycenter)*dsize/dres, dsize/dres),ffdata[0][cutcrop[0]:cutcrop[1],xcenter],ffdata[1][cutcrop[0]:cutcrop[1],xcenter]]

ffcut = ffxcut
x0, x1, y0, y1 = cutcrop[0], cutcrop[1], ycenter, ycenter

sycenter = int(slen/2)+1
sxcenter = int(slen/2)+1

scutcrop = [sxcenter-int(np.ceil(scutnum/2)),sxcenter+int(np.floor(scutnum/2))]

xs = np.add(np.multiply(np.arange(-ssize/2,ssize/2,sres),1),0)

nrightcut = [xs,nrightnv[sycenter,:]]
nleftcut = [xs,nleftnv[sycenter,:]]
blochcut = [xs,blochnv[sycenter,:]]
nlbcut = [xs,nlbnv[sycenter,:]]
nrbcut = [xs,nrbnv[sycenter,:]]

sx0, sx1, sy0, sy1 = scutcrop[0], scutcrop[1], sycenter, sycenter

# ---------------------- BZ LINECUTS -------------------------------------------------
# ------------------------------------------------------------------------------------

bx0, by0 = xcenter, ycenter
simbx0, simby0 = (slen/2), (slen/2)
phinum = 8
bzlclen = 15
dstep = dsize/(dres*sres)
simbzlclen = bzlclen*dstep
bzs = (1e3)*np.arange(-bzlclen*dstep/2, bzlclen*dstep/2-0.000001, dstep)
simbzs = (1e3)*np.arange(-bzlclen*dstep/2, bzlclen*dstep/2, dstep)
bzphi = np.zeros((phinum, bzlclen))
bzphiError = np.zeros((phinum, bzlclen))
simbznlbphi = np.zeros((phinum, bzlclen))
simbznrbphi = np.zeros((phinum, bzlclen))
simbznrightphi = np.zeros((phinum, bzlclen))
simbznleftphi = np.zeros((phinum, bzlclen))
simbzblochphi = np.zeros((phinum, bzlclen))
bzrad = np.zeros((2*phinum, 8))
nlbrad = np.zeros((2*phinum, 8))
nrbrad = np.zeros((2*phinum, 8))
blochrad = np.zeros((2*phinum, 8))
nleftrad = np.zeros((2*phinum, 8))
nrightrad = np.zeros((2*phinum, 8))
radcutdiff = np.zeros((5, 2*phinum))

for i in range(0,phinum):
   phic = i*pi/phinum
   bx1, by1 = bx0-bzlclen*np.cos(phic), by0-bzlclen*np.sin(phic)
   bx2, by2 = bx0+bzlclen*np.cos(phic), by0+bzlclen*np.sin(phic)
   bx, by = np.linspace(bx1, bx2, bzlclen), np.linspace(by1, by2, bzlclen)
   bzphi[i] = ndimage.map_coordinates(np.transpose(bzdata), np.vstack((bx,by)), order=1)
   bzphiError[i] = ndimage.map_coordinates(np.transpose(bzdataError), np.vstack((bx,by)), order=1)
   simbx1, simby1 = simbx0-simbzlclen*np.cos(phic), simby0-simbzlclen*np.sin(phic)
   simbx2, simby2 = simbx0+simbzlclen*np.cos(phic), simby0+simbzlclen*np.sin(phic)
   simbx, simby = np.linspace(simbx1, simbx2, bzlclen), np.linspace(simby1, simby2, bzlclen)
   simbznlbphi[i] = ndimage.map_coordinates(np.transpose(nlb[2]), np.vstack((simbx,simby)), order=1)
   simbznrbphi[i] = ndimage.map_coordinates(np.transpose(nrb[2]), np.vstack((simbx,simby)), order=1)
   simbznrightphi[i] = ndimage.map_coordinates(np.transpose(nright[2]), np.vstack((simbx,simby)), order=1)
   simbznleftphi[i] = ndimage.map_coordinates(np.transpose(nleft[2]), np.vstack((simbx,simby)), order=1)
   simbzblochphi[i] = ndimage.map_coordinates(np.transpose(bloch[2]), np.vstack((simbx,simby)), order=1)
   bzrad[i] = bzphi[i][0:8]
   bzrad[phinum+i] = bzphi[i][7:15]
   nlbrad[i] = simbznlbphi[i][0:8]
   nlbrad[phinum+i] = simbznlbphi[i][7:15]
   nrbrad[i] = simbznrbphi[i][0:8]
   nrbrad[phinum+i] = simbznrbphi[i][7:15]
   blochrad[i] = simbzblochphi[i][0:8]
   blochrad[phinum+i] = simbzblochphi[i][7:15]
   nleftrad[i] = simbznleftphi[i][0:8]
   nleftrad[phinum+i] = simbznleftphi[i][7:15]
   nrightrad[i] = simbznrightphi[i][0:8]
   nrightrad[phinum+i] = simbznrightphi[i][7:15]

for i in range(0, 2*phinum):
    radcutdiff[0][i] = np.sum((nlbrad[i]-bzrad[i])**2)
    radcutdiff[1][i] = np.sum((blochrad[i]-bzrad[i])**2)
    radcutdiff[2][i] = np.sum((nleftrad[i]-bzrad[i])**2)
    radcutdiff[3][i] = np.sum((nrbrad[i]-bzrad[i])**2)
    radcutdiff[4][i] = np.sum((nrightrad[i]-bzrad[i])**2)

# ---------------------- PLOTS -------------------------------------------------------
# ------------------------------------------------------------------------------------

savepath = '/Users/alec/UCSB/papers/irmn/figures/'
my_dpi = 96
size = 800
size_inches = size/my_dpi

plt.close('all')

fig, ax = plt.subplots()
plt.plot(radcutdiff[0], color='#2D7DD2', linewidth=2.0, label=r'$\gamma_h$ = 103.2 degrees')
plt.plot(radcutdiff[1], color='#E33553', linewidth=2.0, label=u'Bloch')
plt.plot(radcutdiff[2], color='#E58A25', linewidth=2.0, label=u'left-handed Néel')
plt.plot(radcutdiff[3], color='#000000', linewidth=2.0, label=r'$\gamma_h$ = 76.8 degrees')
plt.plot(radcutdiff[4], color='#97CC04', linewidth=2.0, label=u'right-handed Néel')
plt.legend(loc=2, borderaxespad=0., prop={'size':12})
plt.ylabel('residual')
plt.xlabel(r'radial cut angle ($\times\pi/8$)')
fp.format_plot(plt, 600, 600, 0, 50)
pylab.savefig(savepath+'residuals.png', format='png')


fig1, ax1 = plt.subplots()
fig1.set_size_inches(6, 6)
ax1.set_axis_off()
im1 = plt.imshow(nleftnv, cmap='bone',interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
# fp.format_plot(plt, 450, 450, 0, 50)
pylab.savefig(savepath+'nv_sim_'+str(filenum)+filespec+'.png', format='png', dpi=my_dpi)


fig1, ax1 = plt.subplots()
fig1.set_size_inches(8, 8)
ax1.set_axis_off()
im1 = plt.imshow(ffdata[0], cmap='bone',interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
pylab.savefig(savepath+'nv_'+str(filenum)+filespec+'.pdf', format='pdf', dpi=my_dpi)


fig1, ax1 = plt.subplots()
fig1.set_size_inches(8, 8)
ax1.set_axis_off()
im1 = plt.imshow(bzdata, cmap='bone',interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
pylab.savefig(savepath+'bz_'+str(filenum)+filespec+'.pdf', format='pdf', dpi=my_dpi)


fig1, ax1 = plt.subplots()
fig1.set_size_inches(6, 6)
ax1.set_axis_off()
im1 = plt.imshow(nlb[2], cmap='bone',interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
pylab.savefig(savepath+'bz_sim_'+str(filenum)+filespec+'.png', format='png', dpi=my_dpi)


fig1, ax1 = plt.subplots()
fig1.set_size_inches(6, 4)
# plt.plot(blochcut[0],blochcut[1],linewidth=2.0, color='#F97304', label='Bloch')
plt.plot(nlbcut[0],nlbcut[1],linewidth=2.0, color='#2D7DD2', label=r'$\gamma_h$ = 157.3 degrees')
# plt.plot(nleftcut[0],nleftcut[1],linewidth=2.0, color='#97CC04', label=u'left-handed Néel')
plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="data")
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
plt.ylim([-5,40])
plt.xlim([-1.0,1.0])
# fp.format_plot(plt, 600, 450, 0, 450)
pylab.savefig(savepath+'linecut_'+str(filenum)+filespec+'.pdf')

fig, axes = plt.subplots(ncols=2, nrows=int(phinum/2), sharex=True, sharey=True)
fig.set_size_inches(8, 4.8)
# xlocator=MaxNLocator(prune='both', nbins=4)
# ylocator=MaxNLocator(prune='both', nbins=6)
for j in range(0,2):
    for i in range(0,int(phinum/2)):
        axes[i,j].plot(simbzs, simbznrightphi[int(i+(phinum/2)*j)], color='#E58A25', linewidth=2.0, label=u'right-handed Néel', zorder=-8)
        axes[i,j].plot(simbzs, simbzblochphi[int(i+(phinum/2)*j)], color='#E33553', linewidth=2.0, label='Bloch', zorder=-7)
        axes[i,j].plot(simbzs, simbznlbphi[int(i+(phinum/2)*j)], color='#2183BC', linewidth=2.0, label=r'$\gamma_h$ = 157.3 degrees', zorder=-5)
        axes[i,j].plot(simbzs, simbznrbphi[int(i+(phinum/2)*j)], color='#000000', linewidth=2.0, label=r'$\gamma_h$ = 22.7 degrees', zorder=-5)
        axes[i,j].plot(simbzs, simbznleftphi[int(i+(phinum/2)*j)], color='#A1C844', linewidth=2.0, label=u'left-handed Néel', zorder=-6)
        (_, caps, _) = axes[i,j].errorbar(bzs, bzphi[int(i+(phinum/2)*j)], bzphiError[int(i+(phinum/2)*j)],color='#000000',linewidth=1.0,fmt='.',label=r"reconstructed B$_z$", zorder=-4)
        for cap in caps:
             cap.set_markeredgewidth(1)
        axes[i,j].tick_params(width=1)
        # axes[i,j].yaxis.set_major_locator(ylocator)
        # axes[i,j].xaxis.set_major_locator(xlocator)
        axes[i,j].set_xticks(np.linspace(-250,250,3))
        axes[i,j].set_yticks(np.linspace(-30,0,3))
pylab.ylim([-40,10])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':12})
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.7, top=0.95, wspace=0, hspace=0)
# pylab.savefig(savepath+'bzcuts_'+str(filenum)+filespec+'.pdf', format='pdf')


plt.show()
