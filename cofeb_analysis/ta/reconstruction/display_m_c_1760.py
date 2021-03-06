# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T16:32:53-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-03-14T00:05:13-05:00



import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy import signal
import format_plot as fp

# mpl.rcParams['toolbar'] = 'None'

pi = np.pi
scannum = 1903
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
my = -np.multiply(np.sin(phi),mr) #negative to compensate for origin at image top

vec_space = 5
x = np.arange(0,len(mz[0]),vec_space)
y = x.copy()
mx_vec = mx[::vec_space,::vec_space].copy()
my_vec = my[::vec_space,::vec_space].copy()

plt.close('all')

fig1, ax1 = plt.subplots(figsize=(3,3))
fig1.set_size_inches(4, 4)
im1 = plt.imshow(mz, cmap='jet')
Q = plt.quiver(x, y, mx_vec, my_vec, units='width', scale=20, width=0.005)
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 290, 290, 0, 0, no_axes = True)
pylab.savefig('/Users/alec/UCSB/scan_images/mz_recon_'+str(scannum)+'.png', format='png')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mr, cmap='jet')
Q = plt.quiver(x, y, mx_vec, my_vec, units='width')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 500, 500, 0, 320, no_axes = True)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(phi_seed, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 290, 290, 290, 0, no_axes = True)

# fig1, ax1 = plt.subplots()
# im1 = plt.imshow(mx_seed, cmap='jet')
# fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
# fp.format_plot(plt, 290, 290, 290, 320, no_axes = True)

# fig1, ax1 = plt.subplots()
# im1 = plt.imshow(my_seed, cmap='jet')
# fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
# fp.format_plot(plt, 290, 290, 290, 660, no_axes = True)

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(phi, cmap='jet', interpolation='Nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 500, 500, 580, 0, no_axes = True)
pylab.savefig('/Users/alec/UCSB/scan_images/phi_recon_'+str(scannum)+'.png', format='png')

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(mx, cmap='jet', interpolation='Nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 500, 500, 400, 400, no_axes = True)
pylab.savefig('/Users/alec/UCSB/scan_images/mx_recon_'+str(scannum)+'.png', format='png')


fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(my, cmap='jet', interpolation='Nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 500, 500, 900, 400, no_axes = True)
pylab.savefig('/Users/alec/UCSB/scan_images/my_recon_'+str(scannum)+'.png', format='png')


plt.show()
