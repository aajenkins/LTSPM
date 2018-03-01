# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T11:41:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-25T11:32:46-06:00



import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy import ndimage

import stray_field_calc_thick as sfct
import calc_NV_field as cNVf
import plotting.format_plots_tkagg as fp

scanSize = 3e-6
DWWidth = 35e-9
radius = 400e-9
Ms = 654577.0
t = 1e-09
helicity = 48*np.pi/180
x = np.linspace(-(scanSize/2), (scanSize/2), 400)
xx, yy = np.meshgrid(x,x)
phi_grid = np.arctan2(yy,xx) + np.pi

mz = np.tanh((np.sqrt(xx**2 + yy**2) - radius)/DWWidth)
mx = np.sqrt(1-mz**2) * np.cos(phi_grid+helicity)
my = np.sqrt(1-mz**2) * np.sin(phi_grid+helicity)

b, scd, vcd, meff = sfct.stray_field_calc_thick(mx, my, mz, Ms, t, 3e-6, 60e-9)

theta = 0.96753783354062239
phi = 0
bNV = cNVf.calc_NV_field(b[0], b[1], b[2]+9.5e-4, theta, phi)

plt.close('all')

fig, ax = plt.subplots()
plt.imshow(mx)
fig, ax = plt.subplots()
plt.imshow(my)
fig, ax = plt.subplots()
plt.imshow(mz)
fig, ax = plt.subplots()
plt.imshow(bNV)
plt.colorbar()
fp.format_plots(plt)

bNVc = np.exp(-(bNV**2)/(2*(2e-4)**2))
fig, ax = plt.subplots()
plt.imshow(bNVc, cmap='bone_r', interpolation='bicubic')
circle1 = plt.Circle((200, 200), (0.4e-6)*len(x)/scanSize, fill=False, color='r')
ax.add_artist(circle1)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)


plt.show()
