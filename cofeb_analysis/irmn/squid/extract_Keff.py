# @Author: Jenkins Alec <alec>
# @Date:   2017-01-23T16:23:24-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-03-28T22:28:19-07:00

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy import stats
from scipy import interpolate

import format_plot as fp

#--------------- SAMPLE PARAMETERS ---------------------------------

sampleAreaCM = 0.26
thicknessCM = 0.911e-7

#--------------- HYSTERESIS FIT ------------------------------------

dataPath = "/Users/alec/UCSB/SQUID_data/irmn-1-02-7-17/"

sq_data_all_ll = np.genfromtxt(dataPath+'IrMn_parallell_Ms4.rso.dat', delimiter=',', skip_header=31)
sq_data_all_norm = np.genfromtxt(dataPath+'IrMn_normal_Ms3.rso.dat', delimiter=',', skip_header=31)

length_ll = 25
sq_data_ll = np.transpose(sq_data_all_ll[0:length_ll, [2, 4]])
length_norm = len(sq_data_all_norm)
sq_data_norm = np.transpose(sq_data_all_norm[0:length_norm, [2, 4]])

baselen_ll = 10
baseline_data_ll = sq_data_ll[[0,1],length_ll-baselen_ll:length_ll]
baselen_norm = 8
baseline_data_norm = sq_data_norm[[0,1],length_norm-baselen_norm:length_norm]

slope_ll, intercept_ll, r_value_ll, p_value_ll, std_err_ll = stats.linregress(baseline_data_ll[0], baseline_data_ll[1])
slope_norm, intercept_norm, r_value_norm, p_value_norm, std_err_norm = stats.linregress(baseline_data_norm[0], baseline_data_norm[1])

sq_data_base_ll = sq_data_ll.copy()
for i in range(0, len(sq_data_ll[0])):
	sq_data_base_ll[1, i] = sq_data_ll[1, i] - (slope_ll*sq_data_ll[0, i])

sq_data_base_norm = sq_data_norm.copy()
for i in range(0, len(sq_data_norm[0])):
	sq_data_base_norm[1, i] = sq_data_norm[1, i] - (slope_norm*sq_data_norm[0, i])

sq_data_base_norm_extend = sq_data_base_norm.copy()
asym_norm = np.mean(sq_data_base_norm[1,length_norm-baselen_norm:length_norm])
asym_ll = np.mean(sq_data_base_ll[1,length_ll-baselen_ll:length_ll])
ms_avg = (asym_norm + asym_ll)/2

asym_norm_std = np.std(sq_data_base_norm[1,length_norm-baselen_norm:length_norm])
asym_ll_std = np.std(sq_data_base_ll[1,length_ll-baselen_ll:length_ll])
ms_avg_std = np.sqrt( asym_norm_std**2 + asym_ll_std**2 )



Mst_avg = ms_avg/sampleAreaCM

print('avg. ms = '+str(ms_avg))
print('Mst_avg = '+str(Mst_avg))
print('ms_avg_std = '+str(ms_avg_std))

sq_data_base_ll[1, :] = np.multiply(sq_data_base_ll[1, :], ms_avg/asym_ll)
sq_data_base_norm_extend[1, :] = np.multiply(sq_data_base_norm_extend[1, :], ms_avg/asym_norm)
sq_data_base_norm_extend = np.transpose(np.concatenate(
	(np.transpose(sq_data_base_norm_extend), np.array([[sq_data_base_ll[0,-1], ms_avg]])),axis=0))

normint = interpolate.splrep(sq_data_base_norm_extend[0], sq_data_base_norm_extend[1], s=0, k = 1)
yint = interpolate.splev(sq_data_base_ll[0],normint)
sq_data_base_norm_int = sq_data_base_ll[0], yint

area_ll = np.trapz(sq_data_base_ll[1], sq_data_base_ll[0])
area_norm = np.trapz(sq_data_base_norm_extend[1], sq_data_base_norm_extend[0])

area_diff = area_norm-area_ll

Keff = area_diff/(sampleAreaCM*thicknessCM)

print('area between norm and ll = '+str(area_diff))
print('Keff = '+str(Keff))


# plt.close('all')

# fig1, ax1 = plt.subplots()
# plt.fill_between(sq_data_base_ll[0], sq_data_base_ll[1], sq_data_base_norm_int[1])
# plt.plot(sq_data_base_norm_int[0], sq_data_base_norm_int[1], 'g', label="normal")
# plt.plot(sq_data_base_ll[0], sq_data_base_ll[1], 'r', label="parallel")
# plt.legend(bbox_to_anchor=(0.9, 0.5), loc=1, borderaxespad=0., prop={'size':12})
# ax1.set_xlabel('H (Oe)')
# ax1.set_ylabel('m (emu)')
# fp.format_plot(plt, 600, 400, 0, 50)
# pylab.savefig('/Users/alec/UCSB/scan_images/irmn/SQUID_measurement.png')

fig1, ax1 = plt.subplots(figsize=(5,4.5))
plt.plot(sq_data_base_norm_extend[0], sq_data_base_norm_extend[1]*(1e6), color="#F97304", label=r"H$_\perp$", linewidth=2.0)
plt.plot(sq_data_base_ll[0], sq_data_base_ll[1]*(1e6), color="#2D7DD2", label="H$_{||}$", linewidth=2.0)
plt.legend(bbox_to_anchor=(0.92, 0.3), loc=1, borderaxespad=0., prop={'size':15})
# ax1.set_xlabel('H (Oe)')
# ax1.set_ylabel('magnetic moment (emu x10^-6)')
ax1.set_xticks(np.linspace(0.0,4000,5))
plt.axis([0,4000,0,30])
# fp.format_plot(plt, 500, 400, 0, 50)
pylab.savefig('/Users/alec/UCSB/papers/irmn/figures/Keff.pdf')

fig1, ax1 = plt.subplots()
plt.plot(sq_data_base_norm[0,], sq_data_base_norm[1])
fp.format_plot(plt, 600, 400, 650, 50)

fig1, ax1 = plt.subplots()
plt.plot(sq_data_ll[0], sq_data_ll[1])
plt.plot(sq_data_norm[0], sq_data_norm[1])
fp.format_plot(plt, 600, 400, 0, 450)

plt.show()
