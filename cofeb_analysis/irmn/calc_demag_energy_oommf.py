# @Author: Jenkins Alec <alec>
# @Date:   2017-02-17T14:20:47-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-21T15:32:31-07:00



import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import imread
from scipy import interpolate
from skimage import measure
import json

import format_plot as fp

pi = np.pi
mu0 = 4*pi*(1e-7)
deadlayer = False

def get_interpolated_contour(contour):
    slen = len(contour[:, 1])
    cint, u = interpolate.splprep([contour[:, 1], contour[:, 0]], s=3)
    xint, yint = interpolate.splev(np.linspace(0,1,5*slen),cint)
    contour_length = 0
    for i in range(0, 5*slen-1):
        contour_length += ( np.sqrt((xint[i+1] - xint[i])**2
                                                  + (yint[i+1] - yint[i])**2) )
    return contour_length, [yint, xint]


# pass parameters with SI units
def get_demag_energy(magn, h, Ms, t, res):

    mlen = len(magn)
    demag_energy = 0
    volume = t*(res**2)

    for j in range(dlen):
        for i in range(dlen):
            demag_energy += magn[j][i] * h[j][i]

    demag_energy = -(1/2) * mu0 * (Ms**2) * volume * demag_energy
    return demag_energy


scannum = 1836
impath = '/Users/alec/UCSB/scan_images/irmn_contour_selection/domains'+str(scannum)+'.png'
domains = imread(impath, flatten=True)
domains = np.add(np.multiply(2/255,domains),-1)

path = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
material_params_path = path+'material_parameters.json'
scan_params_path = path+'scan_parameters.json'
with open(material_params_path, 'r') as fread:
    material_params = json.load(fread)
with open(scan_params_path, 'r') as fread:
    scan_params = json.load(fread)

scansize = 1.0e-5
slen = len(domains)
res = scansize/slen
res_difference = 2e-9

if (deadlayer):
    Ms = material_params['Ms2']
    thickness = material_params['t2']
else:
    Ms = material_params['Ms']
    thickness = material_params['t']

Keff = material_params['Keff']

dlen = len(domains)
dres = 2.0e-9

contours = measure.find_contours(domains, 0.0)
num_contours = len(contours)

total_wall_length = 0

smoothed_contours = np.zeros_like(contours)
contour_length = 0
for i in range(0, len(contours)):
    contour_length, smoothed_contour = get_interpolated_contour(contours[i])
    smoothed_contours[i] = np.transpose(smoothed_contour)
    total_wall_length += res*contour_length

print('total_wall_length = ' + str(total_wall_length))

total_wall_length_norm = total_wall_length/( (res*(dlen-2))**2 )

ohfPath = '/Users/alec/UCSB/oommf/data_and_runs/irmn/'
filename = 'irmn_magn'+str(1+deadlayer)+'-hdemag.ohf'
ohfDemagH = np.genfromtxt(ohfPath+filename)
filenameDres = 'irmn_magn_dres'+str(1+deadlayer)+'-hdemag.ohf'
ohfDemagHDres = np.genfromtxt(ohfPath+filenameDres)

demagH = np.zeros((dlen,dlen))
demagHDres = np.zeros((dlen,dlen))

for j in range(dlen):
    for i in range(dlen):
        demagH[j][i] = ohfDemagH[i+dlen*j][2]
        demagHDres[j][i] = ohfDemagHDres[i+dlen*j][2]


demag_energy0 = get_demag_energy(domains, demagH, Ms, thickness, res)
demag_energy1 = get_demag_energy(domains, demagHDres, Ms, thickness, res+dres)

demag_energy0_norm = demag_energy0/((res*dlen)**2)
demag_energy1_norm = demag_energy1/(((res+dres)*dlen)**2)
delta_demag_energy_norm = (demag_energy1_norm - demag_energy0_norm)/dres
delta_total_wall_length_norm = -total_wall_length_norm/res

domain_wall_energy_density = -( delta_demag_energy_norm /
                               (delta_total_wall_length_norm * thickness) )

print('domain_wall_energy_density = ' + str(domain_wall_energy_density))

As = (domain_wall_energy_density**2)/(16*Keff)
DW_width = domain_wall_energy_density/(4*Keff)

print('DW_width (no pi) = ' + str(DW_width))
print('As = ' + str(As))

plt.close('all')

fig, ax = plt.subplots()
im = plt.imshow(demagH, cmap='jet', interpolation='nearest')
plt.colorbar(im, fraction=0.046, pad=0.04)
fp.format_plot(plt, 500, 500, 50, 50)

fig, ax = plt.subplots()
im = plt.imshow(domains, cmap='jet', interpolation='nearest')
plt.colorbar(im, fraction=0.046, pad=0.04)
fp.format_plot(plt, 500, 500, 550, 50)

plt.show()
