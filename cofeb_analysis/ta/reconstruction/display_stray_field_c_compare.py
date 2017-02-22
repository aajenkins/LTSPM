# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T16:32:53-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-16T12:43:29-08:00



import numpy as np
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc as sfc

pi = np.pi

path = '/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/'
mx = np.loadtxt(path+"mbx_test.dat")
my = np.loadtxt(path+"mby_test.dat")
mz = np.loadtxt(path+"mbz_test.dat")

sfc_result = sfc.stray_field_calc(mx, my, mz, 1.0e-3, 4000, 80e-9)

hz = sfc_result[4][2]
hzc = np.loadtxt(path+'stray_field_test_c.txt', delimiter=',')

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(hz, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(hzc, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

plt.show()
