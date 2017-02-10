# @Author: Jenkins Alec <alec>
# @Date:   2017-01-23T16:23:24-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-08T15:30:55-08:00

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

import format_plot as fp

#--------------- SAMPLE PARAMETERS ---------------------------------

sampleAreaCM = 0.26
thicknessCM = 0.9e-7

#--------------- HYSTERESIS FIT ------------------------------------

datapath = "/Users/alec/UCSB/SQUID_data/irmn-1-02-7-17/"

sq_data_all = np.genfromtxt(datapath+'IrMn_normal_Ms6.rso.dat', delimiter=',', skip_header=31)

length1 = 20
length2 = 41
length3 = 41
sq_data1 = np.transpose(sq_data_all[0:length1, [2, 4]])
sq_data2 = np.transpose(sq_data_all[length1:length1+length2, [2, 4]])
sq_data3 = np.transpose(sq_data_all[length1+length2:length1+length2+length3, [2, 4]])

sq_data = np.multiply(np.add(sq_data2, np.fliplr(sq_data3)), 1/2.)

dlen = length2

baselen = 6
baseline_data1 = sq_data[[0,1],0:baselen]
baseline_data2 = sq_data[[0,1],dlen-baselen:dlen]
#
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(baseline_data1[0], baseline_data1[1])
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(baseline_data2[0], baseline_data2[1])
slope = (slope1 + slope2)/2

sq_data_base = sq_data.copy()
for i in range(0, len(sq_data[0])):
	sq_data_base[1, i] = sq_data[1, i] - (slope*sq_data[0, i])

asym_data1 = sq_data_base[1,0:baselen]
asym_data2 = np.multiply(sq_data_base[1,dlen-baselen:dlen],-1)
asym_data = np.concatenate((asym_data1, asym_data2))

ms = [np.mean(asym_data), np.std(asym_data)]

Mst = ms[0]/sampleAreaCM

Ms = ms[0]/(sampleAreaCM*thicknessCM)

print('ms = '+str(ms[0])+' +/- '+str(ms[1])+' emu')
print('Ms*t (emu/cm^2) = '+str(Mst))
print('Ms*t (A) = '+str(Mst*(1e3)*(1e-2)))
print('Ms (emu/cc) = '+str(Ms))

plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(sq_data1[0], sq_data1[1])
fp.format_plot(plt, 600, 400, 0, 50)

fig1, ax1 = plt.subplots()
plt.plot(sq_data[0], sq_data[1])
fp.format_plot(plt, 600, 400, 650, 50)

fig1, ax1 = plt.subplots()
plt.plot(sq_data_base[0], sq_data_base[1])
fp.format_plot(plt, 600, 400, 0, 450)

plt.show()
