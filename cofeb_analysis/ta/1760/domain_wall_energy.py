# @Author: Jenkins Alec <alec>
# @Date:   2017-01-28T18:37:33-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-14T20:03:06-08:00

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
from scipy.ndimage import imread
from scipy import interpolate
from skimage import measure
from skimage import morphology

import calc_demag_tensor_newell as cdtn
import demag_energy_calc_newell as decn
import format_plot as fp

def get_interpolated_contour(contour):
    slen = len(contour[:, 1])
    cint, u = interpolate.splprep([contour[:, 1], contour[:, 0]], s=3)
    xint, yint = interpolate.splev(np.linspace(0,1,5*slen),cint)
    contour_length = 0
    for i in range(0, 5*slen-1):
        contour_length += np.sqrt((xint[i+1] - xint[i])**2
                                                  + (yint[i+1] - yint[i])**2)
    return contour_length, [yint, xint]

scannum = 809
impath = '/Users/alec/UCSB/cofeb_analysis_data/ta/domains'+str(scannum)+'.png'
domains = imread(impath, flatten=True)
domains = np.add(np.multiply(2/255,domains),-1)

dataPath = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
filespec = 'Msfixed'
cal_params = np.loadtxt(dataPath+'cal_parameters_'+filespec+'.txt', delimiter=',')

demag_tensor = np.loadtxt('/Users/alec/UCSB/cofeb_analysis_data/ta/demag_tensor.txt', delimiter=',')

Ms = cal_params[0]
thickness = cal_params[1]
Keff = cal_params[6]

MsSI = Ms*(1.0e3)
thicknessSI = thickness*(1.0e-2)
KeffSI = Keff*(1.0e-1)

scanvsize = 4
scansize = (5.0e-6)*scanvsize
slen = len(domains)
res = scansize/slen
res_difference = 2.0e-9

contours = measure.find_contours(domains[1:slen-1,1:slen-1], 0.0)
num_contours = len(contours)

total_wall_length = 0

smoothed_contours = np.zeros_like(contours)
contour_length = 0
for i in range(0, len(contours)):
    contour_length, smoothed_contour = get_interpolated_contour(contours[i])
    smoothed_contours[i] = np.transpose(smoothed_contour)
    total_wall_length += res*contour_length

print('total_wall_length = ' + str(total_wall_length))

total_wall_length_norm = total_wall_length/(res*(slen-2))

demag_tensor0 = cdtn.calc_demag_tensor_newell(slen, res, thickness)
demag_energy0, h0 = decn.demag_energy_calc_newell(domains, demag_tensor, MsSI, res, thicknessSI)
demag_energy0_norm = demag_energy0 / ( (res*(slen-2))**2 )
demag_tensor1 = cdtn.calc_demag_tensor_newell(slen, res+res_difference, thickness)
demag_energy1, h1 = decn.demag_energy_calc_newell(domains, demag_tensor, MsSI, res+res_difference, thicknessSI)
demag_energy1_norm = demag_energy1 / ( ((res+res_difference)*(slen-2))**2 )

delta_demag_energy_norm = (demag_energy1_norm - demag_energy0_norm)/res_difference
delta_total_wall_length_norm = -total_wall_length_norm/res

domain_wall_energy_density = -( delta_demag_energy_norm /
                               (delta_total_wall_length_norm * thicknessSI) )

print('domain_wall_energy_density = ' + str(domain_wall_energy_density))

As = (domain_wall_energy_density**2)/(16*KeffSI)
DW_width = domain_wall_energy_density/(4*KeffSI)

print('DW_width (no pi) = ' + str(DW_width))
print('As = ' + str(As))

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(domains[1:slen-1,1:slen-1], cmap='Greys', interpolation='nearest',alpha=0.5)

for n, contour in enumerate(smoothed_contours):
    ax1.plot(contour[:, 1], contour[:, 0], linewidth=1)
ax1.axis('image')
ax1.set_xticks([])
ax1.set_yticks([])
fp.format_plot(plt, 800, 800, 50, 50)
# pylab.savefig('/Users/alec/UCSB/scan_images/irmn/DW_length'+str(scannum)+'.png')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(domains, cmap='Greys', interpolation='nearest')
fp.format_plot(plt, 400, 400, 850, 50)

plt.show()
