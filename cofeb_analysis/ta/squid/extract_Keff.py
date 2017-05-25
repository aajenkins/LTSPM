# @Author: Jenkins Alec <alec>
# @Date:   2017-01-23T16:23:24-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-03-13T16:38:20-05:00

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

dataPath = "/Users/alec/UCSB/SQUID_data/ta-7/thinned/"

sq_data_all_ll = np.genfromtxt(dataPath+'Ta_ll_Ms_thinned_3.rso.dat', delimiter=',', skip_header=31)
sq_data_all_norm = np.genfromtxt(dataPath+'Ta_normal_Ms_thinned.rso.dat', delimiter=',', skip_header=31)
sq_data_all_norm_lowH = np.genfromtxt(dataPath+'Ta_normal_Ms_thinned_lowfield.rso.dat', delimiter=',', skip_header=31)


length_ll = 25
sq_data_ll = np.transpose(sq_data_all_ll[0:length_ll, [2, 4]])
length_norm = len(sq_data_all_norm)-1
sq_data_norm = np.transpose(sq_data_all_norm[0:length_norm, [2, 4]])
length_norm_lowH = len(sq_data_all_norm_lowH)-2
sq_data_norm_lowH = np.transpose(sq_data_all_norm_lowH[0:length_norm_lowH, [2, 4]])

baselen_ll = 10
baseline_data_ll = sq_data_ll[[0,1],length_ll-baselen_ll:length_ll]
baselen_norm = 15
baseline_data_norm = sq_data_norm[[0,1],length_norm-baselen_norm:length_norm]
baselen_norm_lowH = 12
baseline_data_norm_lowH = sq_data_norm_lowH[[0,1],length_norm_lowH-baselen_norm_lowH:length_norm_lowH]

slope_ll, intercept_ll, r_value_ll, p_value_ll, std_err_ll = stats.linregress(baseline_data_ll[0], baseline_data_ll[1])
slope_norm, intercept_norm, r_value_norm, p_value_norm, std_err_norm = stats.linregress(baseline_data_norm[0], baseline_data_norm[1])
slope_norm_lowH, intercept_norm_lowH, r_value_norm_lowH, p_value_norm_lowH, std_err_norm_lowH = stats.linregress(baseline_data_norm_lowH[0], baseline_data_norm_lowH[1])


sq_data_base_ll = sq_data_ll.copy()
for i in range(0, length_ll):
	sq_data_base_ll[1, i] = sq_data_ll[1, i] - (slope_ll*sq_data_ll[0, i])
# sq_data_base_ll[1][0] = 0.0

sq_data_base_norm = sq_data_norm.copy()
for i in range(0, length_norm):
	sq_data_base_norm[1, i] = sq_data_norm[1, i] - (slope_norm*sq_data_norm[0, i])
# sq_data_base_norm[1][0] = 0.0

sq_data_base_norm_lowH = sq_data_norm_lowH.copy()
for i in range(0, length_norm_lowH):
	sq_data_base_norm_lowH[1, i] = sq_data_norm_lowH[1, i] - (slope_norm_lowH*sq_data_norm_lowH[0, i])


# sq_data_base_norm_extend = sq_data_base_norm.copy()
asym_norm_lowH = np.mean(sq_data_base_norm_lowH[1,length_norm_lowH-baselen_norm_lowH:length_norm_lowH])
asym_norm = np.mean(sq_data_base_norm[1,length_norm-baselen_norm:length_norm])
asym_ll = np.mean(sq_data_base_ll[1,length_ll-baselen_ll:length_ll])
ms_avg = (asym_norm + asym_ll)/2

# mst std
asym_norm_lowH_std = np.std(sq_data_base_norm_lowH[1,length_norm_lowH-baselen_norm_lowH:length_norm_lowH])

Mst_avg = ms_avg/sampleAreaCM
Mst_normal = asym_norm/sampleAreaCM
Mst_normal_lowH = asym_norm_lowH/sampleAreaCM
Mst_normal_lowH_std = asym_norm_lowH_std/sampleAreaCM

sq_data_base_norm[1] = sq_data_base_norm[1] * asym_norm_lowH/asym_norm
sq_data_base_ll[1] = sq_data_base_ll[1] * asym_norm_lowH/asym_ll
sq_data_base_norm_extend = np.zeros((2,length_norm+1))
sq_data_base_norm_extend[0] = np.append(sq_data_base_norm[0], sq_data_base_ll[0,-1])
sq_data_base_norm_extend[1] = np.append(sq_data_base_norm[1], sq_data_base_ll[1,-1])

print('avg. ms = '+str(ms_avg))
print('Mst_avg = '+str(Mst_avg))
print('Mst_normal = '+str(Mst_normal))
print('Mst_normal_lowH = '+str(Mst_normal_lowH))
print('Mst_normal_lowH_std = '+str(Mst_normal_lowH_std))


area_ll = np.trapz(sq_data_base_ll[1], sq_data_base_ll[0])
area_norm = np.trapz(sq_data_base_norm_extend[1], sq_data_base_norm_extend[0])

area_diff = area_norm-area_ll

Keff = area_diff/(sampleAreaCM*thicknessCM)

print('area between norm and ll = '+str(area_diff))
print('Keff = '+str(Keff))


plt.close('all')

fig1, ax1 = plt.subplots(figsize=(5,4.5))
plt.plot(sq_data_base_norm_extend[0], sq_data_base_norm_extend[1]*(1e6), color="#F97304", label=r"H$_\perp$", linewidth=2.0)
plt.plot(sq_data_base_ll[0], sq_data_base_ll[1]*(1e6), color="#2D7DD2", label="H$_{||}$", linewidth=2.0)
plt.legend(bbox_to_anchor=(0.92, 0.3), loc=1, borderaxespad=0., prop={'size':15})
# ax1.set_xlabel('H (Oe)')
# ax1.set_ylabel('magnetic moment (emu x10^-6)')
ax1.set_yticks(np.linspace(0.0,4.0,5))
plt.axis([0,1000,0,4])
fp.format_plot(plt, 500, 400, 0, 50)
pylab.savefig('/Users/alec/UCSB/SQUID_data/ta-7/Keff_thinned.pdf')

fig1, ax1 = plt.subplots()
plt.plot(sq_data_base_norm[0], sq_data_base_norm[1])
fp.format_plot(plt, 500, 400, 500, 50)

fig1, ax1 = plt.subplots()
plt.plot(sq_data_base_norm_lowH[0,], sq_data_base_norm_lowH[1])
fp.format_plot(plt, 500, 400, 1000, 50)

fig1, ax1 = plt.subplots()
plt.plot(sq_data_ll[0], sq_data_ll[1])
plt.plot(sq_data_norm[0], sq_data_norm[1])
fp.format_plot(plt, 600, 400, 0, 450)

plt.show()
