# @Author: Jenkins Alec <alec>
# @Date:   2017-01-20T18:13:22-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-16T15:13:26-08:00



import numpy as np
import matplotlib.pyplot as plt
import stray_field_calc as sfc
import format_plot as fp

heightnm = 40
height = heightnm*(1.0e-7)
Ms = 6.496e2
t = 1.0e-7
Mst = Ms*t
fsize = 4.0e-6
MstSInm = Ms*(1e3)

baspath = "/Users/alec/UCSB/cofeb_analysis_data/stray_field_test/"
filenames = ["mnrx_test","mnry_test","mnrz_test","mnlx_test","mnly_test","mnlz_test"]
numfiles = len(filenames)
m = []
for i in range(0,numfiles):
    m.append(np.transpose(np.loadtxt(baspath+filenames[i]+".dat")))

sdr, vcdr, meffr, hkr, hr = sfc.stray_field_calc(m[1],m[0],m[2],MstSInm,2500,heightnm)
sdl, vcdl, meffl, hkl, hl = sfc.stray_field_calc(m[4],m[3],m[5],MstSInm,2500,heightnm)

hintx = np.loadtxt("/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/neelright_x_"+str(heightnm)+"_30._400_10..dat")
hinty = np.loadtxt("/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/neelright_y_"+str(heightnm)+"_30._400_10..dat")
hintz = np.loadtxt("/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/neelright_z_"+str(heightnm)+"_30._400_10..dat")
hintx = np.multiply(hintx, (1e-3)*MstSInm)
hinty = np.multiply(hinty, (1e-3)*MstSInm)
hintz = np.multiply(hintz, (1e-3)*MstSInm)

flen = len(hr[0])
ilen = len(hintx)

ires = 10.0
fres = 2500/flen

ix = np.arange(-ilen*ires/2, ilen*ires/2, ires)
fx = np.arange(-flen*fres/2, flen*fres/2, fres)

plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(ix, hintx, 'r-')
plt.plot(fx, hr[0][int(flen/2), :], 'b-')

fig1, ax1 = plt.subplots()
plt.plot(ix, hintz, 'r-')
plt.plot(fx, hr[2][int(flen/2), :], 'b-')
fp.format_plot(plt, 650, 400, 0, 450)

plt.show()
