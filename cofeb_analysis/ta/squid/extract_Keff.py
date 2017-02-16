# @Author: Jenkins Alec <alec>
# @Date:   2017-01-23T16:23:24-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-14T21:46:50-08:00

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy import stats
from scipy import interpolate

import format_plot as fp

#--------------- SAMPLE PARAMETERS ---------------------------------

sampleAreaCM = 0.055
thicknessCM = 1.0e-7

#--------------- HYSTERESIS FIT ------------------------------------

dataPath = "/Users/alec/UCSB/SQUID_data/ta-7/"

sq_data_all_ll = np.genfromtxt(dataPath+'Ta_parallell_Ms_back_subtract.rso.dat', delimiter=',', skip_header=31)
sq_data_all_norm = np.genfromtxt(dataPath+'Ta_normal_Ms2.rso.dat', delimiter=',', skip_header=31)

length_ll = 8
sq_data_ll = np.transpose(sq_data_all_ll[0:length_ll, [2, 4]])
length_norm = len(sq_data_all_norm)
sq_data_norm = np.transpose(sq_data_all_norm[0:length_norm, [2, 4]])

baselen_ll = 3
baseline_data_ll = sq_data_ll[[0,1],length_ll-baselen_ll:length_ll]
baselen_norm = 12
baseline_data_norm = sq_data_norm[[0,1],length_norm-baselen_norm:length_norm]

slope_ll, intercept_ll, r_value_ll, p_value_ll, std_err_ll = stats.linregress(baseline_data_ll[0], baseline_data_ll[1])
slope_norm, intercept_norm, r_value_norm, p_value_norm, std_err_norm = stats.linregress(baseline_data_norm[0], baseline_data_norm[1])

sq_data_base_ll = sq_data_ll.copy()
for i in range(0, len(sq_data_ll[0])):
	sq_data_base_ll[1, i] = sq_data_ll[1, i] - (slope_ll*sq_data_ll[0, i])
sq_data_base_ll[1][0] = 0.0

sq_data_base_norm = sq_data_norm.copy()
for i in range(0, len(sq_data_norm[0])):
	sq_data_base_norm[1, i] = sq_data_norm[1, i] - (slope_norm*sq_data_norm[0, i])
sq_data_base_norm[1][0] = 0.0

sq_data_base_norm_extend = sq_data_base_norm.copy()
asym_norm = np.mean(sq_data_base_norm[1,length_norm-baselen_norm:length_norm])
asym_ll = np.mean(sq_data_base_ll[1,length_ll-baselen_ll:length_ll])
ms_avg = (asym_norm + asym_ll)/2

Mst_avg = ms_avg/sampleAreaCM

print('avg. ms = '+str(ms_avg))
print('Mst_avg = '+str(Mst_avg))

sq_data_base_ll[1, :] = np.multiply(sq_data_base_ll[1, :], ms_avg/asym_ll)
sq_data_base_norm_extend[1, :] = np.multiply(sq_data_base_norm_extend[1, :], ms_avg/asym_norm)
sq_data_base_norm_extend = np.transpose(np.concatenate(
	(np.transpose(sq_data_base_norm_extend), np.array([[sq_data_base_ll[0,-1], ms_avg]])),axis=0))

# normint = interpolate.splrep(sq_data_base_norm_extend[0], sq_data_base_norm_extend[1], s=0, k = 1)
# yint = interpolate.splev(sq_data_base_ll[0],normint)
# sq_data_base_norm_int = sq_data_base_ll[0], yint

area_ll = np.trapz(sq_data_base_ll[1], sq_data_base_ll[0])
area_norm = np.trapz(sq_data_base_norm_extend[1], sq_data_base_norm_extend[0])

area_diff = area_norm-area_ll

Keff = area_diff/(sampleAreaCM*thicknessCM)

print('area between norm and ll = '+str(area_diff))
print('Keff = '+str(Keff))


plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(sq_data_base_norm_extend[0], sq_data_base_norm_extend[1], 'g', label="normal")
plt.plot(sq_data_base_ll[0], sq_data_base_ll[1], 'r', label="parallel")
plt.legend(bbox_to_anchor=(0.9, 0.5), loc=1, borderaxespad=0., prop={'size':12})
ax1.set_xlabel('H (Oe)')
ax1.set_ylabel('m (emu)')
fp.format_plot(plt, 600, 400, 0, 50)
# pylab.savefig('/Users/alec/UCSB/scan_images/irmn/SQUID_measurement.png')

fig1, ax1 = plt.subplots()
plt.plot(sq_data_base_norm[0,], sq_data_base_norm[1])
fp.format_plot(plt, 600, 400, 650, 50)

fig1, ax1 = plt.subplots()
plt.plot(sq_data_ll[0], sq_data_ll[1])
plt.plot(sq_data_norm[0], sq_data_norm[1])
fp.format_plot(plt, 600, 400, 0, 450)

plt.show()
