import numpy as np
import matplotlib.pyplot as plt
import stray_field_calc as sfc
import format_plot as fp

height = 40
Ms = 6.496e5
t = 1
Mst = Ms*t

baspath = "/Users/alec/UCSB/LTSPM/cofeb_analysis/ta/stray_field_sim/"
filenames = ["mbx","mby","mbz"]
numfiles = len(filenames)
m = []
for i in range(0,numfiles):
    m.append(np.transpose(np.loadtxt(baspath+filenames[i]+".dat")))

scd, vcd, meff, hk, h = sfc.stray_field_calc(m[0],m[1],m[2],Mst,2500,height)

hintx = np.loadtxt("/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/bloch_x_"+str(height)+"_30._400_10..dat")
hinty = np.loadtxt("/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/bloch_y_"+str(height)+"_30._400_10..dat")
hintz = np.loadtxt("/Users/alec/UCSB/mathematica/CoFeB-MgO/linecut_simulations/linecut_30_h/bloch_z_"+str(height)+"_30._400_10..dat")
hintx = np.multiply(hintx, (1e-3)*Mst)
hinty = np.multiply(hinty, (1e-3)*Mst)
hintz = np.multiply(hintz, (1e-3)*Mst)

flen = len(h[0])
ilen = len(hintx)

ires = 10
fres = 2500/flen

ix = np.arange(-ilen*ires/2, ilen*ires/2, ires)
fx = np.arange(-flen*fres/2, flen*fres/2, fres)

plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(ix, hintx, 'r-')
plt.plot(fx, h[0][int(flen/2), :], 'b-')
fp.format_plot(plt, 650, 400, 0, 50)

fig1, ax1 = plt.subplots()
plt.plot(ix, hintz, 'r-')
plt.plot(fx, h[2][int(flen/2), :], 'b-')
fp.format_plot(plt, 650, 400, 0, 450)

plt.show()