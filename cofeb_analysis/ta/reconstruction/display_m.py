# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T16:32:53-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-10T15:35:11-08:00



import numpy as np
import matplotlib.pyplot as plt
import format_plot as fp

pi = np.pi

path = '/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/'
mz = np.loadtxt(path+'mzdata_test.txt', delimiter=',')
phi = np.loadtxt(path+'phi_test.txt', delimiter=',')
mr = np.sqrt(1-mz**2)
mx = np.multiply(np.cos(phi),mr)
my = np.multiply(np.sin(phi),mr)

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mz, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mr, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 450, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(phi, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 450)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mx, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 450, 450)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(my, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 850, 450)


plt.show()
