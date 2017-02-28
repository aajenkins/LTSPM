# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:17:32 2016

@author: alec
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
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

thetas = [theta-thetaError,theta,theta+thetaError]

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

bloch = [[],[],[]]
nleft = [[],[],[]]
nright = [[],[],[]]

for i in range(0, len(errnames)):
    blochfms[i] = np.append(blochfms[i],glob.glob(simpath+'b_*'+errnames[i]+'*_lowres*'+filespec+'.txt'))
    nleftfms[i] = np.append(nleftfms[i],glob.glob(simpath+'nl_*'+errnames[i]+'*lowres*'+filespec+'.txt'))
    nrightfms[i] = np.append(nrightfms[i],glob.glob(simpath+'nr_*'+errnames[i]+'*lowres*'+filespec+'.txt'))
    bloch[i] = [[],[],[]]
    nleft[i] = [[],[],[]]
    nright[i] = [[],[],[]]
    for k in range(0, 3):
        bloch[i][k] = np.loadtxt(blochfms[i][k], delimiter=',')
    bloch[i][2] = bloch[i][2]+zfields
    for k in range(0, 3):
        nleft[i][k] = np.loadtxt(nleftfms[i][k], delimiter=',')
    nleft[i][2] = nleft[i][2]+zfields
    for k in range(0, 3):
        nright[i][k] = np.loadtxt(nrightfms[i][k], delimiter=',')
    nright[i][2] = nright[i][2]+zfields

#simulation constants
ssize = 2.5
slen = len(bloch[0][0][0])
sres = ssize/slen

blochnv = np.zeros((3,slen,slen))
nleftnv = np.zeros((3,slen,slen))
nrightnv = np.zeros((3,slen,slen))


for j in range(0, slen):
    for i in range(0, slen):
        blochnv[1][j,i] = cNV.calc_NV_field(bloch[1][0][j,i], bloch[1][1][j,i], bloch[1][2][j,i], thetas[1], phi)
        nleftnv[1][j,i] = cNV.calc_NV_field(nleft[1][0][j,i], nleft[1][1][j,i], nleft[1][2][j,i], thetas[1], phi)
        nrightnv[1][j,i] = cNV.calc_NV_field(nright[1][0][j,i], nright[1][1][j,i], nright[1][2][j,i], thetas[1], phi)

        blochnv[0][j,i] = max(cNV.calc_NV_field(bloch[0][0][j,i], bloch[0][1][j,i], bloch[0][2][j,i], thetas[0], phi),
                                     cNV.calc_NV_field(bloch[0][0][j,i], bloch[0][1][j,i], bloch[0][2][j,i], thetas[1], phi),
                                     cNV.calc_NV_field(bloch[0][0][j,i], bloch[0][1][j,i], bloch[0][2][j,i], thetas[2], phi))
        nleftnv[0][j,i] = max(cNV.calc_NV_field(nleft[0][0][j,i], nleft[0][1][j,i], nleft[0][2][j,i], thetas[0], phi),
                                     cNV.calc_NV_field(nleft[0][0][j,i], nleft[0][1][j,i], nleft[0][2][j,i], thetas[1], phi),
                                     cNV.calc_NV_field(nleft[0][0][j,i], nleft[0][1][j,i], nleft[0][2][j,i], thetas[2], phi))
        nrightnv[0][j,i] = max(cNV.calc_NV_field(nright[0][0][j,i], nright[0][1][j,i], nright[0][2][j,i], thetas[0], phi),
                                      cNV.calc_NV_field(nright[0][0][j,i], nright[0][1][j,i], nright[0][2][j,i], thetas[1], phi),
                                      cNV.calc_NV_field(nright[0][0][j,i], nright[0][1][j,i], nright[0][2][j,i], thetas[2], phi))

        blochnv[2][j,i] = min(cNV.calc_NV_field(bloch[2][0][j,i], bloch[2][1][j,i], bloch[2][2][j,i], thetas[0], phi),
                                     cNV.calc_NV_field(bloch[2][0][j,i], bloch[2][1][j,i], bloch[2][2][j,i], thetas[1], phi),
                                     cNV.calc_NV_field(bloch[2][0][j,i], bloch[2][1][j,i], bloch[2][2][j,i], thetas[2], phi))
        nleftnv[2][j,i] = min(cNV.calc_NV_field(nleft[2][0][j,i], nleft[2][1][j,i], nleft[2][2][j,i], thetas[0], phi),
                                     cNV.calc_NV_field(nleft[2][0][j,i], nleft[2][1][j,i], nleft[2][2][j,i], thetas[1], phi),
                                     cNV.calc_NV_field(nleft[2][0][j,i], nleft[2][1][j,i], nleft[2][2][j,i], thetas[2], phi))
        nrightnv[2][j,i] = min(cNV.calc_NV_field(nright[2][0][j,i], nright[2][1][j,i], nright[2][2][j,i], thetas[0], phi),
                                      cNV.calc_NV_field(nright[2][0][j,i], nright[2][1][j,i], nright[2][2][j,i], thetas[1], phi),
                                      cNV.calc_NV_field(nright[2][0][j,i], nright[2][1][j,i], nright[2][2][j,i], thetas[2], phi))

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,maxfgrad=20)
# misc.imsave('/Users/alec/UCSB/scan_images/full-field/ff'+str(filenum)+'.png', ffdata[0])
bxdata = np.loadtxt(datapath+'bx_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bydata = np.loadtxt(datapath+'by_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')
bzdata = np.loadtxt(datapath+'bz_'+str(filenum)+'_'+filespec+'.txt', delimiter=',')

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

xs = np.add(np.multiply(np.arange(-ssize/2,ssize/2,sres),1),0)

for k in range(0,3):
    bcut[k] = [xs,blochnv[k][sycenter, :]]
    nrcut[k] = [xs,nrightnv[k][sycenter, :]]
    nlcut[k] = [xs,nleftnv[k][sycenter, :]]
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
bzlclen = 20
simbzlclen = int(bzlclen*dsize/(dres*sres))
bzphi = np.zeros((phinum, bzlclen))
simbzphi = np.zeros((3, phinum, simbzlclen))
for i in range(0,phinum):
   phic = i*pi/phinum
   bx1, by1 = bx0-bzlclen*np.cos(phic), by0-bzlclen*np.sin(phic)
   bx2, by2 = bx0+bzlclen*np.cos(phic), by0+bzlclen*np.sin(phic)
   bx, by = np.linspace(bx1, bx2, bzlclen), np.linspace(by1, by2, bzlclen)
   bzs = np.arange(0,bzlclen*dsize/dres,dsize/dres)
   bzphi[i] = ndimage.map_coordinates(np.transpose(bzdata), np.vstack((bx,by)), order=1)

   simbx1, simby1 = simbx0-simbzlclen*np.cos(phic), simby0-simbzlclen*np.sin(phic)
   simbx2, simby2 = simbx0+simbzlclen*np.cos(phic), simby0+simbzlclen*np.sin(phic)
   simbx, simby = np.linspace(simbx1, simbx2, simbzlclen), np.linspace(simby1, simby2, simbzlclen)
   simbzs = np.arange(0,simbzlclen*sres,sres)
   simbzphi[0][i] = ndimage.map_coordinates(np.transpose(nright[1][2]), np.vstack((simbx,simby)), order=1)
   simbzphi[1][i] = ndimage.map_coordinates(np.transpose(bloch[1][2]), np.vstack((simbx,simby)), order=1)
   simbzphi[2][i] = ndimage.map_coordinates(np.transpose(nleft[1][2]), np.vstack((simbx,simby)), order=1)

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
# for i in range(0,phinum):
#    phic = i*pi/phinum
#    bx1, by1 = bx0-bzlclen*np.cos(phic), by0-bzlclen*np.sin(phic)
#    bx2, by2 = bx0+bzlclen*np.cos(phic), by0+bzlclen*np.sin(phic)
# plt.plot([bx1, bx2], [by1, by2], 'r-')
pylab.savefig('/Users/alec/UCSB/scan_images/bzcut_diagram_'+str(filenum)+filespec+'.tif', format='tif', dpi=my_dpi)

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
# plt.subplots_adjust(left=0.05, bottom=0.05, right=0.7, top=0.95, wspace=0, hspace=0)
# pylab.savefig('/Users/alec/UCSB/scan_images/bzcuts_'+str(filenum)+filespec+'.pdf', format='pdf')
#
#
fig1, ax1 = plt.subplots()
plt.fill_between(bcut[0][0], bcut[0][1], bcut[2][1],color='#96bee8',linewidth=1.0)
plt.fill_between(nlcut[0][0], nlcut[0][1], nlcut[2][1],color='#fcb981',linewidth=1.0)
plt.fill_between(nrcut[0][0], nrcut[0][1], nrcut[2][1],color='#cbe581',linewidth=1.0)
plt.plot(bcut[1][0],bcut[1][1],color='#2D7DD2',linewidth=2.0,label="Bloch")
plt.plot(nlcut[1][0],nlcut[1][1],color='#F97304',linewidth=2.0, label=u'left-handed Néel')
plt.plot(nrcut[1][0],nrcut[1][1],color='#97CC04',linewidth=2.0, label=u'right-handed Néel')
plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="NV ESR data")
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
pylab.ylim([0,40])
pylab.xlim([-1.2,1.2])
fp.format_plot(plt, 600, 450, 0, 450)
pylab.savefig('/Users/alec/UCSB/scan_images/linecut_'+str(filenum)+filespec+'.pdf', format='pdf')

plt.show()
