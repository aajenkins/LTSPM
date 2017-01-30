import numpy as np
import matplotlib.pyplot as plt
import stray_field_calc as sfc
import format_plot as fp

scannum = 1760
path = '/Users/alec/UCSB/LTSPM/cofeb_analysis/ta/'+str(scannum)+'/'
filespec = 'Msnotfixed'
cal_params = np.loadtxt(path+'cal_parameters_'+filespec+'.txt', delimiter=',')

Ms = cal_params[0]
t = cal_params[1]
theta = cal_params[2]
phi = cal_params[3]
height = cal_params[4]

Mst = Ms*t

baspath = "/Users/alec/UCSB/LTSPM/cofeb_analysis/ta/1760/stray_field_sim/"
filenames = ["mbx_Msnotfixed","mby_Msnotfixed","mbz_Msnotfixed"]
filenames_symmetric = ["mbx_symmetric","mby_symmetric","mbz_symmetric"]
numfiles = len(filenames)
m = []
for i in range(0,numfiles):
    m.append(np.loadtxt(baspath+filenames[i]+".dat"))

scd, vcd, meff, hk, h = sfc.stray_field_calc(m[0],m[1],m[2],Mst,2500,height)

m_symmetric = []
for i in range(0,numfiles):
    m_symmetric.append(np.loadtxt(baspath+filenames_symmetric[i]+".dat"))

scd, vcd, meff, hk, h_symmetric = sfc.stray_field_calc(m_symmetric[0],m_symmetric[1],m_symmetric[2],Mst,2500,height)

flen = len(h[0])
fres = 2500/flen
fx = np.arange(-flen*fres/2, flen*fres/2, fres)

plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(fx, m_symmetric[0][int(flen/2), :], 'r-')
plt.plot(fx, m[0][int(flen/2), :], 'b-')
fp.format_plot(plt, 650, 400, 0, 50)

fig1, ax1 = plt.subplots()
plt.plot(fx, h_symmetric[2][int(flen/2), :], 'r-')
plt.plot(fx, h[2][int(flen/2), :], 'b-')
fp.format_plot(plt, 650, 400, 0, 450)

fig1, ax1 = plt.subplots()
plt.plot(fx, h_symmetric[2][:, int(flen/2)], 'r-')
plt.plot(fx, h[2][:, int(flen/2)], 'b-')
fp.format_plot(plt, 650, 400, 0, 450)

plt.show()
