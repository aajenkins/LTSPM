# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T16:32:53-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-03-13T14:56:49-05:00



import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy import signal
import format_plot as fp
import fourier_image as fi

# mpl.rcParams['toolbar'] = 'None'

pi = np.pi
scannum = 1903
mgauge = 'neel_left'
path = '/Users/alec/UCSB/mathematica/magnetization_reconstruction/'
mz = np.loadtxt(path+'mz_1760_'+mgauge+'_gauge.dat')[1:-1,1:-1]
phi = np.loadtxt(path+'phi_1760_'+mgauge+'_gauge.dat')
phiSeed = np.loadtxt(path+'phiSeed_1760_'+mgauge+'_gauge.dat')
mr = np.sqrt(1-mz**2)
mz = fi.window_image(mz, power=1/4)
mr = fi.window_image(mr, power=1/4)
mx = np.multiply(np.cos(phi),mr)
my = -np.multiply(np.sin(phi),mr) #negative to compensate for origin at image top

vec_space = 2
x = np.arange(0,len(mz[0]),vec_space)
y = x.copy()
mx_vec = mx[::vec_space,::vec_space].copy()
my_vec = my[::vec_space,::vec_space].copy()

plt.close('all')

fig1, ax1 = plt.subplots(figsize=(3,3))
fig1.set_size_inches(4, 4)
im1 = plt.imshow(mz, cmap='jet', interpolation='Nearest')
Q = plt.quiver(x, y, mx_vec, my_vec, units='width', scale=20, width=0.005)
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 500, 500, 0, 0, no_axes = True)
# pylab.savefig('/Users/alec/UCSB/scan_images/mz_neel_'+mgauge+'_'+str(scannum)+'.png', format='png')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(mr, cmap='jet', interpolation='Nearest')
Q = plt.quiver(x, y, mx_vec, my_vec, units='width', scale=20, width=0.005)
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 500, 500, 0, 320, no_axes = True)

fig, ax = plt.subplots()
plt.plot(mz[28,:])
fp.format_plot(plt, 600, 400, 290, 320)
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

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(phiSeed, cmap='jet', interpolation='Nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 500, 500, 580, 500, no_axes = True)
# pylab.savefig('/Users/alec/UCSB/scan_images/phi_'+mgauge+'_'+str(scannum)+'.png', format='png')

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(mx, cmap='jet', interpolation='Nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 500, 500, 400, 400, no_axes = True)
# pylab.savefig('/Users/alec/UCSB/scan_images/mx_'+mgauge+'_'+str(scannum)+'.png', format='png')


fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 4)
im1 = plt.imshow(my, cmap='jet', interpolation='Nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 500, 500, 900, 400, no_axes = True)
# pylab.savefig('/Users/alec/UCSB/scan_images/my_'+mgauge+'_'+str(scannum)+'.png', format='png')


plt.show()
