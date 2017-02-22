# @Author: Jenkins Alec <alec>
# @Date:   2017-01-20T18:13:22-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-18T14:41:45-08:00



import numpy as np
import matplotlib.pyplot as plt
import stray_field_calc_thick as sfct
import format_plot as fp

import fourier_image as fi

heightnm = 40
height = heightnm*(1.0e-7)
Ms = 6.496e2
t = 1.0e-7
Mst = Ms*t
fsize = 4.0e-6
MstSInm = Ms*(1e3)

def edge_field(x, *args):
    h = args[0]
    x0 = 0
    bx = (2e-3)*MstSInm*(h/(h**2+(x-x0)**2))
    bz = -(2e-3)*MstSInm*((x-x0)/(h**2+(x-x0)**2))

    return bx, bz

baspath = "/Users/alec/UCSB/cofeb_analysis_data/stray_field_test/"
filenames = ["medgex_test","medgey_test","medgez_test"]
numfiles = len(filenames)
m = []
for i in range(0,numfiles):
    m.append(np.transpose(np.loadtxt(baspath+filenames[i]+".dat")))

mwx = fi.window_image(m[1])
mwy = fi.window_image(m[0])
mwz = fi.window_image(m[2])
sdr, vcdr, meffr, hkr, h = sfct.stray_field_calc_thick(mwx,mwy,mwz,Mst,2.5e-4,height)

flen = len(h[0])

fres = 2500/flen
fx = np.arange(-flen*fres/2, flen*fres/2, fres)


bx, bz = edge_field(fx, heightnm)

plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(fx, h[1][:,int(flen/2)], 'b-')
plt.plot(fx, bx, 'r-')
fp.format_plot(plt, 650, 400, 0, 450)

fig1, ax1 = plt.subplots()
plt.plot(fx, h[2][:,int(flen/2)], 'b-')
plt.plot(fx, bz, 'r-')
fp.format_plot(plt, 650, 400, 0, 450)

plt.show()
