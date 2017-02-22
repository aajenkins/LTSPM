# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T16:32:53-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-15T20:16:44-08:00



import numpy as np
import matplotlib.pyplot as plt
import format_plot as fp

pi = np.pi

path = '/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/'
mz = np.loadtxt(path+'mzdata_test.txt')
phi = np.loadtxt(path+'phi_c_test.txt', delimiter=',')
phi_seed = np.loadtxt(path+'phi_seed_c_test.txt', delimiter=',')
phi_seed = phi_seed%(2*pi)
phi = phi%(2*pi)
mr = np.sqrt(1-mz**2)
mx_seed = np.multiply(np.cos(phi_seed),mr)
my_seed = np.multiply(np.sin(phi_seed),mr)
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
im1 = plt.imshow(phi_seed, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 450, 450)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mx, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 250, 250, 1100, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(my, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 250, 250, 1100, 450)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mx_seed, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 250, 250, 850, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(my_seed, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 250, 250, 850, 450)


plt.show()
