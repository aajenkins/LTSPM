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
sfieldpre = 1
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
simpath = datapath+'stray_field_sim/dw_set/'
errnames = ["lower", "mean", "upper"]
filespec = "Msfixed"

dw0s = np.arange(50,200,20)
numdw = len(dw0s)
nleftfms = []
nleft = []
nlcut = []
for i in range(0, numdw):
    nleftfms.append([[],[],[]])
    nleft.append([[],[],[]])
    nlcut.append([[],[],[]])

for j in range(0, numdw):
    for i in range(0, len(errnames)):
        nleftfms[j][i] = np.append(nleftfms[j][i],glob.glob(simpath+'nl_*'+errnames[i]+'*lowres_1903_'+str(dw0s[j])+'.txt'))
        nleft[j][i] = [[],[],[]]
        for k in range(0, 3):
            nleft[j][i][k] = np.loadtxt(nleftfms[j][i][k], delimiter=',')

#simulation constants
ssize = 2.5
slen = len(nleft[0][0][0][0])
sres = ssize/slen

nleftnv = np.zeros((numdw, 3, slen, slen))

for m in range(0, numdw):
    for k in range(0,3):
        for j in range(0, slen):
            for i in range(0, slen):
                nleftnv[m][k][j,i] = np.multiply(sfieldpre,np.abs((nleft[m][k][0][j,i]*np.sin(theta)*np.cos(phi))+
                                                             (nleft[m][k][1][j,i]*np.sin(theta)*np.sin(phi))+
                                                             (np.add(nleft[m][k][2][j,i],zfields)*np.cos(theta))))

ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(filenum)+'-esrdata/fitdata.txt',dres,dres,maxfgrad=30)

ycenter = 27
xcenter = 22

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

for m in range(0, numdw):
    for k in range(0,3):
        nlcut[m][k] = [xs,nleftnv[m][k][sycenter,:]]

sx0, sx1, sy0, sy1 = scutcrop[0], scutcrop[1], sycenter, sycenter

# ---------------------- PLOTS -------------------------------------------------------
# ------------------------------------------------------------------------------------

plt.close('all')

fig1, ax1 = plt.subplots()
# plt.fill_between(nlcut[dw][0][0], nlcut[dw][0][1], nlcut[dw][2][1],color='#F97304',alpha=0.5,linewidth=1.0)
for dw in range(0,len(dw0s)):
    plt.plot(nlcut[dw][1][0],nlcut[dw][1][1],linewidth=2.0, label=str(dw0s[dw]))
# plt.plot(nlcut[dw][1][0],nlcut[dw][1][1],color='#F97304',linewidth=2.0, label=u'left-handed Néel')
plt.errorbar(ffcut[0],ffcut[1],yerr=ffcut[2],color='#ED1035',fmt='.',label="data")
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
pylab.xlim([-1.0,1.0])
fp.format_plot(plt, 600, 450, 0, 450)
pylab.savefig('/Users/alec/UCSB/scan_images/linecut_'+str(filenum)+filespec+'.png')

plt.show()
