# @Author: Jenkins Alec <alec>
# @Date:   2017-01-28T18:37:33-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-04T20:57:36-08:00


import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import imread
from scipy import interpolate
from skimage import measure

import demag_energy_calc_newell as decn
import format_plot as fp

def get_contour_length(contour):
    slen = len(contour[:, 1])
    cint, u = interpolate.splprep([contour[:, 1], contour[:, 0]], s=0)
    xint, yint = interpolate.splev(np.linspace(0,1,5*slen),cint)
    contour_length = 0
    for i in range(0, 5*slen-1):
        contour_length += np.sqrt((xint[i+1] - xint[i])**2
                                                  + (yint[i+1] - yint[i])**2)
    return contour_length

scannum = 1836
impath = '/Users/alec/UCSB/scan_images/irmn/domains'+str(scannum)+'.png'
domains = imread(impath, flatten=True)
domains = np.add(np.multiply(2/255,domains),-1)

scanvsize = 2
scansize = (5e-6)*scanvsize
slen = len(domains)
res = scansize/slen
res_difference = 1e-9
thickness = 0.911e-9
Ms = 1.044e6

contours = measure.find_contours(domains[1:slen-1,1:slen-1], 0.0)
num_contours = len(contours)

total_wall_length = 0

for n, contour in enumerate(contours):
    total_wall_length += res*get_contour_length(contour)

print('total_wall_length = ' + str(total_wall_length))

total_wall_length_norm = total_wall_length/(res*(slen-2))

demag_energy0, h0, d0 = decn.demag_energy_calc_newell(domains, Ms, res, thickness)
demag_energy0_norm = demag_energy0 / ( (res*(slen-2))**2 )
demag_energy1, h1, d1 = decn.demag_energy_calc_newell(domains, Ms, res+res_difference, thickness)
demag_energy1_norm = demag_energy1 / ( ((res+res_difference)*(slen-2))**2 )

delta_demag_energy_norm = (demag_energy1 - demag_energy0)/res_difference
delta_total_wall_length_norm = -total_wall_length_norm/res

domain_wall_energy_density = -( delta_demag_energy_norm /
                               (delta_total_wall_length_norm * thickness) )

print('domain_wall_energy_density = ' + str(domain_wall_energy_density))

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(domains[1:slen-1,1:slen-1], cmap='Greys', interpolation='nearest',alpha=0.5)

for n, contour in enumerate(contours):
    ax1.plot(contour[:, 1], contour[:, 0], linewidth=1)
ax1.axis('image')
ax1.set_xticks([])
ax1.set_yticks([])

fp.format_plot(plt, 400, 400, 50, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(domains, cmap='Greys', interpolation='nearest')
fp.format_plot(plt, 400, 400, 450, 50)

plt.show()
