# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T16:32:53-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-12T20:12:45-08:00



import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import signal
import format_plot as fp

# mpl.rcParams['toolbar'] = 'None'

pi = np.pi

path = '/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/'
mz = np.loadtxt(path+'mzdata_norm.txt')
phi = np.loadtxt(path+'phi_c_1760.txt', delimiter=',')
phi_seed = np.loadtxt(path+'phi_seed_c_1760.txt', delimiter=',')
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
fp.format_plot(plt, 290, 290, 0, 0, no_axes = True)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mr, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 290, 290, 0, 320, no_axes = True)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(phi_seed, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 290, 290, 290, 0, no_axes = True)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mx_seed, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 290, 290, 290, 320, no_axes = True)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(my_seed, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 290, 290, 290, 660, no_axes = True)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(phi, cmap='jet', interpolation='Nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 500, 500, 580, 0, no_axes = True)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mx, cmap='jet', interpolation='Nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 500, 500, 0, 320, no_axes = True)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(my, cmap='jet', interpolation='Nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 500, 500, 580, 660, no_axes = True)


plt.show()
