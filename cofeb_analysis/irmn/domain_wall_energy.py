# @Author: Jenkins Alec <alec>
# @Date:   2017-01-28T18:37:33-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-01-29T21:09:59-08:00


import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import imread
from scipy import interpolate
from skimage import measure

import demag_energy_calc as dec
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

impath = '/Users/alec/UCSB/scan_images/irmn/domains1835.png'
domains = imread(impath, flatten=True)

scanvsize = 2
scansize = (5e-6)*scanvsize
res = scansize/len(domains)
res_difference = 1e-9
thickness = 0.911e-9
Ms = 1.044e6

domains = np.add(np.multiply(2/255,domains),-1)

contours = measure.find_contours(domains, 0.0)
num_contours = len(contours)

total_wall_length = 0

for n, contour in enumerate(contours):
    total_wall_length += get_contour_length(contour)

print('total_wall_length = ' + str(total_wall_length*res))

demag_energy0 = dec.get_demag_energy(domains, Ms, res, thickness)
demag_energy1 = dec.get_demag_energy(domains, Ms, res+res_difference, thickness)

d_demag_energy = (demag_energy1 - demag_energy0)/res_difference
d_total_wall_length = total_wall_length/res

domain_wall_energy_density = d_demag_energy/(d_total_wall_length * thickness)

print('domain_wall_energy_density = ' + str(domain_wall_energy_density))

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(domains, cmap='Greys', interpolation='nearest',alpha=0.5)

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
