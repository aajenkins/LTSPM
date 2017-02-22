# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T16:32:53-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-16T13:37:42-08:00



import numpy as np
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc as sfc
import stray_field_calc_discrete as sfcd

pi = np.pi

scannum = 1760
path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
filespec = 'Msfixed'
cal_params = np.loadtxt(path+'cal_parameters_'+filespec+'.txt', delimiter=',')

Ms = cal_params[0]
t = cal_params[1]
theta = cal_params[2]
phi = cal_params[3]
height = cal_params[4]

Mst = Ms*t

baspath = "/Users/alec/UCSB/cofeb_analysis_data/ta/stray_field_sim/"
filenames = ["mbx_Msfixed","mby_Msfixed","mbz_Msfixed"]
numfiles = len(filenames)
m = []
for i in range(0,numfiles):
    m.append(np.loadtxt(baspath+filenames[i]+".dat"))

print('calc discrete')
scd, vcd, meff, hk, h = sfcd.stray_field_calc_discrete(m[0],m[1],m[2],Mst,2500,height)
print('calc cont')
scd_d, vcd_d, meff_d, hk_d, h_d = sfc.stray_field_calc(m[0],m[1],m[2],Mst,2500,height)

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(h[2], cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(h_d[2], cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

fig1, ax1 = plt.subplots()
im1 = plt.plot(h[2][250,:])
im1 = plt.plot(h_d[2][250,:])
fp.format_plot(plt, 400, 400, 450, 50)


plt.show()
