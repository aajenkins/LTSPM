# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T11:41:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-09-18T14:54:09-07:00



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
DWWidth2 = 15e-9
radius = 400e-9
Ms = 654577.0
t = 1e-09
helicity = 90*np.pi/180
height = 1e-7
x = np.linspace(-(scanSize/2), (scanSize/2), 400)
xx, yy = np.meshgrid(x,x)
phi_grid = np.arctan2(yy,xx) + np.pi

mz = np.tanh((np.sqrt(xx**2 + yy**2) - radius)/DWWidth)
mx = np.sqrt(1-mz**2) * np.cos(phi_grid+helicity)
my = np.sqrt(1-mz**2) * np.sin(phi_grid+helicity)

mz2 = np.tanh((np.sqrt(xx**2 + yy**2) - radius)/DWWidth2)
mx2 = np.sqrt(1-mz**2) * np.cos(phi_grid+helicity)
my2 = np.sqrt(1-mz**2) * np.sin(phi_grid+helicity)

b, scd, vcd, meff = sfct.stray_field_calc_thick(mx, my, mz, Ms, t, 3e-6, height)
b2, scd2, vcd2, meff2 = sfct.stray_field_calc_thick(mx2, my2, mz2, Ms, t, 3e-6, height)

theta = 0.96753783354062239
phi = np.pi
taField = 7.0e-4
wField = 7.0e-4
coptField = 52.0e-4
width = 10e6

bNVta = cNVf.calc_NV_field(b[0], b[1], b[2]+taField, theta, phi)
bNVw = cNVf.calc_NV_field(b[0], b[1], b[2]+wField, theta, phi)
bNVcopt = cNVf.calc_NV_field(b2[0], b2[1], b2[2]+coptField, theta, phi)

plt.close('all')

bNVtac = np.divide(1, (bNVta)**2 + (width/(2*(28e9)))**2)
bNVwc = np.divide(1, (bNVw)**2 + (width/(2*(28e9)))**2)
bNVcoptc = np.divide(1, (bNVcopt-coptField*np.cos(theta))**2 + (width/(2*(28e9)))**2)

fig, ax = plt.subplots()
plt.imshow(bNVtac, cmap='bone_r', interpolation='bicubic')
circle1 = plt.Circle((200, 200), (0.4e-6)*len(x)/scanSize, fill=False, color='r')
ax.add_artist(circle1)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)

fig, ax = plt.subplots()
plt.imshow(bNVwc, cmap='bone_r', interpolation='bicubic')
circle1 = plt.Circle((200, 200), (0.4e-6)*len(x)/scanSize, fill=False, color='r')
ax.add_artist(circle1)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)

fig, ax = plt.subplots()
plt.imshow(bNVcoptc, cmap='bone_r', interpolation='bicubic')
circle1 = plt.Circle((200, 200), (0.4e-6)*len(x)/scanSize, fill=False, color='r')
ax.add_artist(circle1)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)


plt.show()
