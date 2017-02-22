# @Author: Jenkins Alec <alec>
# @Date:   2017-01-22T16:34:48-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-18T20:10:11-08:00



import numpy as np
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc_thick as sfct

pi = np.pi

def Bz_analytic(x, Ms, l, z):
    return 2*pi*Ms*np.exp(-2*pi*z/l)*np.cos(2*pi*x/l)


path = '/Users/alec/UCSB/cofeb_analysis_data/stray_field_test/mzStripes.dat'

sz = np.loadtxt(path)
sx = np.zeros_like(sz)
sy = np.zeros_like(sz)

scansize = 2.0e-3

height = 5.0e-5

scd, vcd, meff, hk, h = sfct.stray_field_calc_thick(sy, sx, sz, 1.0, scansize,
                                    height, windowData=True, windowPower=1/16)

hlen = len(h[2])

x = np.linspace(0, scansize, hlen, endpoint=False)
bza = Bz_analytic(x, 1.0, 1.0e-4, height)

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(h[2], cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 500, 500, 0, 50)

fig1, ax1 = plt.subplots()
plt.plot(x, bza)
plt.plot(x, h[2][int(hlen/2),:])
fp.format_plot(plt, 500, 500, 0, 50)

fig1, ax1 = plt.subplots()
plt.plot(x, bza)
fp.format_plot(plt, 500, 500, 500, 50)


plt.show()
