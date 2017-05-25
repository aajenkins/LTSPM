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

height = 149
nleft = [[],[],[]]
nlcut = [[],[],[]]

nleftfms = glob.glob(simpath+'nl_*'+str(height)+'*lowres_1903.txt')
for k in range(0, 3):
    nleft[k] = np.loadtxt(nleftfms[k], delimiter=',')

#simulation constants
ssize = 2.5
slen = len(nleft[0][0])
sres = ssize/slen

nleftnv = np.zeros((3, slen, slen))

for j in range(0, slen):
    for i in range(0, slen):
        nleftThetaError = sfieldpre*gbnve.get_bnv_theta_error(nleft[0][j,i], nleft[1][j,i], nleft[2][j,i], theta, thetaError, phi)
        nleftnv[1][j,i] = sfieldpre*cNV.calc_NV_field(nleft[0][j,i], nleft[1][j,i], nleft[2][j,i], theta, phi)
        nleftnv[0][j,i] = max(nleftnv[1][j,i] - nleftThetaError, 0)
        nleftnv[2][j,i] = nleftnv[1][j,i] + nleftThetaError

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,maxfgrad=30)

ycenter = 28
xcenter = 22

ffdatacrop = ffdata[0][ycenter-20:ycenter+21,xcenter-20:xcenter+21]

nleftnvcrop = [nleftnv[i][5:46,5:46] for i in range(0,len(nleftnv))]

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

for m in range(0, 3):
    nlcut[m] = [xs,nleftnv[m][sycenter,:]]

sx0, sx1, sy0, sy1 = scutcrop[0], scutcrop[1], sycenter, sycenter

# ---------------------- PLOTS -------------------------------------------------------
# ------------------------------------------------------------------------------------

plt.close('all')

fig1, ax1 = plt.subplots()
plt.fill_between(nlcut[0][0], nlcut[0][1], nlcut[2][1],color='#F97304',alpha=0.5,linewidth=1.0)
plt.plot(nlcut[1][0],nlcut[1][1],linewidth=2.0, color='#2D7DD2', label=u'left-handed Néel')
plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="data")
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
plt.ylim([-5,40])
plt.xlim([-1.0,1.0])
fp.format_plot(plt, 600, 450, 0, 450)
pylab.savefig('/Users/alec/UCSB/scan_images/linecut_'+str(filenum)+filespec+'.png')



plt.show()
