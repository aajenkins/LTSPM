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
matplotlib.rcParams['mathtext.default'] = 'regular'

#material parameters
pi = np.pi
scannum = 1760

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
filespec = 'Msfixed'
cal_params_path = path+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

phi = cal_params['phi']
theta = cal_params['theta']
thetaError = cal_params['thetaError']

zfields = 9.5

#file constants
hnum = 1
rnum = 1
dres = 50
dsize = 0.6*5
filenum = 1760

datapath = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
simpath = datapath+'stray_field_sim/'
errnames = ["lower", "mean", "upper"]

blochfms = [[],[],[]]
nleftfms = [[],[],[]]
nrightfms = [[],[],[]]
bnrfms = [[],[],[]]

bloch = [[],[],[]]
nleft = [[],[],[]]
nright = [[],[],[]]
bnr = [[],[],[]]

for i in range(0, len(errnames)):
    blochfms[i] = np.append(blochfms[i],glob.glob(simpath+'b_*'+errnames[i]+'*_lowres*'+filespec+'.txt'))
    nleftfms[i] = np.append(nleftfms[i],glob.glob(simpath+'nl_*'+errnames[i]+'*lowres*'+filespec+'.txt'))
    nrightfms[i] = np.append(nrightfms[i],glob.glob(simpath+'nr_*'+errnames[i]+'*lowres*'+filespec+'.txt'))
    bnrfms[i] = np.append(bnrfms[i],glob.glob(simpath+'h_*'+errnames[i]+'*lowres*'+filespec+'.txt'))
    bloch[i] = [[],[],[]]
    nleft[i] = [[],[],[]]
    nright[i] = [[],[],[]]
    bnr[i] = [[],[],[]]
    for k in range(0, 3):
        bloch[i][k] = np.loadtxt(blochfms[i][k], delimiter=',')
    bloch[i][2] = bloch[i][2]+zfields
    for k in range(0, 3):
        nleft[i][k] = np.loadtxt(nleftfms[i][k], delimiter=',')
    nleft[i][2] = nleft[i][2]+zfields
    for k in range(0, 3):
        nright[i][k] = np.loadtxt(nrightfms[i][k], delimiter=',')
    nright[i][2] = nright[i][2]+zfields
    for k in range(0, 3):
        bnr[i][k] = np.loadtxt(bnrfms[i][k], delimiter=',')
    bnr[i][2] = bnr[i][2]+zfields

#simulation constants
ssize = dsize
slen = len(bloch[0][0][0])
sres = ssize/slen

blochnv = np.zeros((3,slen,slen))
nleftnv = np.zeros((3,slen,slen))
nrightnv = np.zeros((3,slen,slen))
bnrnv = np.zeros((3,slen,slen))


for j in range(0, slen):
    for i in range(0, slen):
        blochThetaError = gbnve.get_bnv_theta_error(bloch[1][0][j,i], bloch[1][1][j,i], bloch[1][2][j,i], theta, thetaError, phi)
        nleftThetaError = gbnve.get_bnv_theta_error(nleft[1][0][j,i], nleft[1][1][j,i], nleft[1][2][j,i], theta, thetaError, phi)
        nrightThetaError = gbnve.get_bnv_theta_error(nright[1][0][j,i], nright[1][1][j,i], nright[1][2][j,i], theta, thetaError, phi)
        bnrThetaError = gbnve.get_bnv_theta_error(bnr[1][0][j,i], bnr[1][1][j,i], bnr[1][2][j,i], theta, thetaError, phi)

        blochnv[1][j,i] = cNV.calc_NV_field(bloch[1][0][j,i], bloch[1][1][j,i], bloch[1][2][j,i], theta, phi)
        nleftnv[1][j,i] = cNV.calc_NV_field(nleft[1][0][j,i], nleft[1][1][j,i], nleft[1][2][j,i], theta, phi)
        nrightnv[1][j,i] = cNV.calc_NV_field(nright[1][0][j,i], nright[1][1][j,i], nright[1][2][j,i], theta, phi)
        bnrnv[1][j,i] = cNV.calc_NV_field(bnr[1][0][j,i], bnr[1][1][j,i], bnr[1][2][j,i], theta, phi)

        blochnv[0][j,i] = cNV.calc_NV_field(bloch[0][0][j,i], bloch[0][1][j,i], bloch[0][2][j,i], theta, phi) + blochThetaError
        nleftnv[0][j,i] = cNV.calc_NV_field(nleft[0][0][j,i], nleft[0][1][j,i], nleft[0][2][j,i], theta, phi) + nleftThetaError
        nrightnv[0][j,i] = cNV.calc_NV_field(nright[0][0][j,i], nright[0][1][j,i], nright[0][2][j,i], theta, phi) + nrightThetaError
        bnrnv[0][j,i] = cNV.calc_NV_field(bnr[0][0][j,i], bnr[0][1][j,i], bnr[0][2][j,i], theta, phi) + bnrThetaError

        blochnv[2][j,i] = max(cNV.calc_NV_field(bloch[2][0][j,i], bloch[2][1][j,i], bloch[2][2][j,i], theta, phi) - blochThetaError,0)
        nleftnv[2][j,i] = max(cNV.calc_NV_field(nleft[2][0][j,i], nleft[2][1][j,i], nleft[2][2][j,i], theta, phi) - nleftThetaError,0)
        nrightnv[2][j,i] = max(cNV.calc_NV_field(nright[2][0][j,i], nright[2][1][j,i], nright[2][2][j,i], theta, phi) - nrightThetaError,0)
        bnrnv[2][j,i] = max(cNV.calc_NV_field(bnr[2][0][j,i], bnr[2][1][j,i], bnr[2][2][j,i], theta, phi) - bnrThetaError,0)

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,maxfgrad=20)
# misc.imsave('/Users/alec/UCSB/scan_images/full-field/ff'+str(filenum)+'.png', ffdata[0])
bxdata = np.loadtxt(datapath+'bx_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bydata = np.loadtxt(datapath+'by_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bzdata = np.loadtxt(datapath+'bz_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bzdataError = np.loadtxt(datapath+'bzError_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')

ycenter = 29
xcenter = 24

cutlength = 2.2

dcutnum = dres*cutlength/dsize
scutnum = cutlength/sres

cutcrop = [xcenter-int(np.ceil(dcutnum/2)), xcenter+int(np.floor(dcutnum/2))]

xd = np.arange(-xcenter*dsize/dres, (dres-xcenter)*dsize/dres, dsize/dres)
xdc = np.arange((cutcrop[0]-xcenter)*dsize/dres, (cutcrop[1]-xcenter)*dsize/dres, dsize/dres)
ffxcut = [xdc, ffdata[0][ycenter,cutcrop[0]:cutcrop[1]],ffdata[1][ycenter,cutcrop[0]:cutcrop[1]]]
ffycut = [np.arange((cutcrop[0]-ycenter)*dsize/dres, (cutcrop[1]-ycenter)*dsize/dres, dsize/dres),ffdata[0][cutcrop[0]:cutcrop[1],xcenter],ffdata[1][cutcrop[0]:cutcrop[1],xcenter]]

ffcut = ffxcut
x0, x1, y0, y1 = cutcrop[0], cutcrop[1], ycenter, ycenter

# ffcut = ffycut
# x0, x1, y0, y1 = xcenter, xcenter, 1, dres-1

sycenter = int(slen/2)
sxcenter = int(slen/2)

sxcutcrop = [sxcenter-int(np.ceil(scutnum/2)), sxcenter+int(np.floor(scutnum/2))]
sycutcrop = [sycenter-int(np.ceil(scutnum/2)), sycenter+int(np.floor(scutnum/2))]

bcut = [[],[],[]]
nrcut = [[],[],[]]
nlcut = [[],[],[]]
bnrcut = [[],[],[]]

xs = np.add(np.multiply(np.arange(-ssize/2,ssize/2,sres),1),0)

for k in range(0,3):
    bcut[k] = [xs,blochnv[k][sycenter, :]]
    nrcut[k] = [xs,nrightnv[k][sycenter, :]]
    nlcut[k] = [xs,nleftnv[k][sycenter, :]]
    bnrcut[k] = [xs,bnrnv[k][sycenter, :]]

    # bcut[k] = [xs,blochnv[k][:, sxcenter]]
    # nrcut[k] = [xs,nrightnv[k][:, sxcenter]]
    # nlcut[k] = [xs,nleftnv[k][:, sxcenter]]

sx0, sx1, sy0, sy1 = sxcutcrop[0], sxcutcrop[1], sycenter, sycenter
# sx0, sx1, sy0, sy1 = sxcenter, sxcenter, sycutcrop[0], sycutcrop[1]

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
simbzblochphi = np.zeros((3, phinum, simbzlclen))
simbznleftphi = np.zeros((3, phinum, simbzlclen))
simbznrightphi = np.zeros((3, phinum, simbzlclen))

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
   simbznrightphi[0][i] = ndimage.map_coordinates(np.transpose(nright[0][2]), np.vstack((simbx,simby)), order=1)
   simbznrightphi[1][i] = ndimage.map_coordinates(np.transpose(nright[1][2]), np.vstack((simbx,simby)), order=1)
   simbznrightphi[2][i] = ndimage.map_coordinates(np.transpose(nright[2][2]), np.vstack((simbx,simby)), order=1)
   simbzblochphi[0][i] = ndimage.map_coordinates(np.transpose(bloch[0][2]), np.vstack((simbx,simby)), order=1)
   simbzblochphi[1][i] = ndimage.map_coordinates(np.transpose(bloch[1][2]), np.vstack((simbx,simby)), order=1)
   simbzblochphi[2][i] = ndimage.map_coordinates(np.transpose(bloch[2][2]), np.vstack((simbx,simby)), order=1)
   simbznleftphi[0][i] = ndimage.map_coordinates(np.transpose(nleft[0][2]), np.vstack((simbx,simby)), order=1)
   simbznleftphi[1][i] = ndimage.map_coordinates(np.transpose(nleft[1][2]), np.vstack((simbx,simby)), order=1)
   simbznleftphi[2][i] = ndimage.map_coordinates(np.transpose(nleft[2][2]), np.vstack((simbx,simby)), order=1)

# ---------------------- PLOTS -------------------------------------------------------
# ------------------------------------------------------------------------------------

plt.close('all')

# display plots
# fig1, ax1 = plt.subplots()
# im1 = plt.imshow(blochnv[1], cmap='bone')
# fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
# plt.plot([sx0, sx1], [sy0, sy1], 'r-')
# plt.axis('image')
# ax1.set_axis_off()
# fp.format_plot(plt, 450, 450, 0, 50)
#
# fig1, ax1 = plt.subplots()
# im1 = plt.imshow(ffdata[0], cmap='bone', interpolation='nearest')
# fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
# plt.plot([x0, x1], [y0, y1], 'r-')
# plt.axis('image')
# ax1.set_axis_off()
# fp.format_plot(plt, 450, 450, 450, 50)
#
# fig1, ax1 = plt.subplots()
# im1 = plt.imshow(bzdata, cmap='bone', interpolation='nearest')
# fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
# plt.axis('image')
# fp.format_plot(plt, 450, 450, 850, 50)
# for i in range(0,phinum):
#    phic = i*2*pi/phinum
#    bx1, by1 = bx0-bzlclen*np.cos(phic), by0-bzlclen*np.sin(phic)
#    bx2, by2 = bx0+bzlclen*np.cos(phic), by0+bzlclen*np.sin(phic)
#    plt.plot([bx1, bx2], [by1, by2], 'r-')
# plt.axis('image')
# fig, axes = plt.subplots(ncols=2, nrows=int(phinum/2), sharex=True, sharey=True)
# for j in range(0,2):
#     for i in range(0,int(phinum/2)):
#         axes[i,j].plot(bzs, bzphi[int(i+(phinum/2)*j)],'r.')
#         axes[i,j].plot(simbzs, simbzphi[0][int(i+(phinum/2)*j)], color='#97CC04', linewidth=2.0, label=u'right-handed Néel')
#         axes[i,j].plot(simbzs, simbzphi[1][int(i+(phinum/2)*j)], color='#2D7DD2', linewidth=2.0, label="Bloch")
#         axes[i,j].plot(simbzs, simbzphi[2][int(i+(phinum/2)*j)], color='#F97304', linewidth=2.0, label=u'left-handed Néel')
#         axes[i,j].get_yaxis().set_visible(False)
#         axes[i,j].text(0.05,.5,u'ϕ = '+'{:2.1f}'.format((i+(phinum/2)*j)/phinum)+' π',
#             horizontalalignment='left', verticalalignment='center',
#             transform=axes[i,j].transAxes, fontsize=12)
# pylab.ylim([-35,20])
# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':12})
# plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)
# fp.format_plot(plt, 900, 900, 450, 50, tight=False)
# plt.subplots_adjust(left=0.05, bottom=0.05, right=0.7, top=0.95, wspace=0, hspace=0)
#
# fig1, ax1 = plt.subplots()
# plt.fill_between(bcut[0][0], bcut[0][1], bcut[2][1],color='#2D7DD2',alpha=0.5,linewidth=1.0)
# plt.fill_between(nlcut[0][0], nlcut[0][1], nlcut[2][1],color='#F97304',alpha=0.5,linewidth=1.0)
# plt.fill_between(nrcut[0][0], nrcut[0][1], nrcut[2][1],color='#97CC04',alpha=0.5,linewidth=1.0)
# plt.plot(bcut[1][0],bcut[1][1],color='#2D7DD2',linewidth=2.0,label="Bloch")
# plt.plot(nlcut[1][0],nlcut[1][1],color='#F97304',linewidth=2.0, label=u'left-handed Néel')
# plt.plot(nrcut[1][0],nrcut[1][1],color='#97CC04',linewidth=2.0, label=u'right-handed Néel')
# plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="data")
# # plt.plot(xd, np.abs(np.sin(-theta)*bxdata[ycenter, :]+np.cos(-theta)*np.add(bzdata[ycenter, :],9.5)))
# plt.legend(loc=1,borderaxespad=1,prop={'size':10})
# pylab.ylim([0,40])
# pylab.xlim([-1.2,1.2])
# fp.format_plot(plt, 600, 450, 0, 450)
#
# plt.show()

# save figures
my_dpi = 96
size = 800
size_inches = size/my_dpi
fig = plt.figure(frameon=False)
fig.set_size_inches(size_inches, size_inches)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
plt.pcolormesh(blochnv[1], cmap='bone')
plt.gca().invert_yaxis()
# plt.plot([sx0, sx1], [sy0, sy1], 'r-')
plt.axis('image')
pylab.savefig('/Users/alec/UCSB/scan_images/ff_sim_'+str(filenum)+filespec+'_bloch.tif', format='tif', dpi=my_dpi)

fig = plt.figure(frameon=False)
fig.set_size_inches(size_inches, size_inches)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
plt.pcolormesh(ffdata[0], cmap='bone')
plt.gca().invert_yaxis()
# plt.plot([x0, x1], [y0, y1], 'r-')
pylab.savefig('/Users/alec/UCSB/scan_images/ff_'+str(filenum)+filespec+'.tif', format='tif', dpi=my_dpi)

fig = plt.figure(frameon=False)
fig.set_size_inches(size_inches, size_inches)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
plt.pcolormesh(bzdata, cmap='bone')
plt.gca().invert_yaxis()
for i in range(0,phinum):
   phic = i*pi/phinum
   bx1, by1 = bx0-bzlclen*np.cos(phic), by0-bzlclen*np.sin(phic)
   bx2, by2 = bx0+bzlclen*np.cos(phic), by0+bzlclen*np.sin(phic)
# plt.plot([bx1, bx2], [by1, by2], 'r-', linewidth=2.0)
pylab.savefig('/Users/alec/UCSB/scan_images/bzcut_diagram_'+str(filenum)+filespec+'.tif', format='tif', dpi=my_dpi)

fig, axes = plt.subplots(ncols=2, nrows=int(phinum/2), sharex=True, sharey=True)
fig.set_size_inches(5, 5)
xlocator=MaxNLocator(prune='both', nbins=6)
ylocator=MaxNLocator(prune='both', nbins=6)
for j in range(0,2):
    for i in range(0,int(phinum/2)):
        axes[i,j].fill_between(simbzs, simbznrightphi[0][int(i+(phinum/2)*j)],simbznrightphi[2][int(i+(phinum/2)*j)], color='#97CC04', zorder=-10)
        axes[i,j].plot(simbzs, simbznrightphi[1][int(i+(phinum/2)*j)], color='#97CC04', linewidth=2.0, label=u'right-handed Néel', zorder=-9)
        axes[i,j].fill_between(simbzs, simbzblochphi[0][int(i+(phinum/2)*j)], simbzblochphi[2][int(i+(phinum/2)*j)], color='#2D7DD2', zorder=-8)
        axes[i,j].plot(simbzs, simbzblochphi[1][int(i+(phinum/2)*j)], color='#2D7DD2', linewidth=2.0, label="Bloch", zorder=-7)
        axes[i,j].fill_between(simbzs, simbznleftphi[0][int(i+(phinum/2)*j)], simbznleftphi[2][int(i+(phinum/2)*j)], color='#F97304', zorder=-6)
        axes[i,j].plot(simbzs, simbznleftphi[0][int(i+(phinum/2)*j)], color='#F97304', linewidth=2.0, label=u'left-handed Néel', zorder=-5)
        (_, caps, _) = axes[i,j].errorbar(bzs, bzphi[int(i+(phinum/2)*j)], bzphiError[int(i+(phinum/2)*j)],color='#ED1035',linewidth=1.0,fmt='.',label=r"reconstructed B$_z$", zorder=-4)
        for cap in caps:
             cap.set_markeredgewidth(1)
        axes[i,j].yaxis.set_major_locator(ylocator)
        axes[i,j].xaxis.set_major_locator(xlocator)
        # if (i != 0):
        #     axes[i,j].get_yaxis().set_visible(False)
        # axes[i,j].text(0.05,.5,u'ϕ = '+'{:2.1f}'.format((i+(phinum/2)*j)/phinum)+' π',
        #     horizontalalignment='left', verticalalignment='center',
        #     transform=axes[i,j].transAxes, fontsize=12)
pylab.ylim([-35,20])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':12})
# plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.95, wspace=0, hspace=0)
pylab.savefig('/Users/alec/UCSB/scan_images/bzcuts_'+str(filenum)+filespec+'.pdf', format='pdf')
pylab.savefig('/Users/alec/UCSB/scan_images/bzcuts_'+str(filenum)+filespec+'.tif', format='tif', dpi=my_dpi)


fig1, ax1 = plt.subplots()
fig1.set_size_inches(6, 5)
plt.fill_between(nlcut[0][0], nlcut[0][1], nlcut[2][1],color='#F97304',linewidth=1.0, zorder=-10)
plt.plot(nlcut[1][0],nlcut[1][1],color='#F97304',linewidth=2.0, label=u'left-handed Néel', zorder=-9)
plt.fill_between(bcut[0][0], bcut[0][1], bcut[2][1],color='#2D7DD2',linewidth=1.0, zorder=-8)
plt.plot(bcut[1][0],bcut[1][1],color='#2D7DD2',linewidth=2.0,label="Bloch", zorder=-7)
plt.fill_between(nrcut[0][0], nrcut[0][1], nrcut[2][1],color='#97CC04',linewidth=1.0, zorder=-6)
plt.plot(nrcut[1][0],nrcut[1][1],color='#97CC04',linewidth=2.0, label=u'right-handed Néel', zorder=-5)
# plt.fill_between(bnrcut[0][0], bnrcut[0][1], bnrcut[2][1],color='#000000',linewidth=1.0, zorder=-6)
# plt.plot(bnrcut[1][0],bnrcut[1][1],color='#000000',linewidth=2.0, label=u'mixture', zorder=-5)
(_, caps, _) = plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="NV ESR data",linewidth=1.0, zorder=-4)
for cap in caps:
     cap.set_markeredgewidth(1)
plt.legend(loc=1,borderaxespad=1,prop={'size':12})
pylab.ylim([0,40])
pylab.xlim([-1.2,1.2])
ax1.set_yticks(np.linspace(0.0,40.0,5))
fp.format_plot(plt, 600, 450, 0, 450)
# pylab.savefig('/Users/alec/UCSB/scan_images/linecut_'+str(filenum)+filespec+'.pdf', format='pdf')
# pylab.savefig('/Users/alec/UCSB/scan_images/linecut_'+str(filenum)+filespec+'.tif', format='tif', dpi=my_dpi)

plt.show()
