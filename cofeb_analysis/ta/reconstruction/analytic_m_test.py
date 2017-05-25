# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T11:41:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-03-14T00:00:59-05:00



import numpy as np
import json
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc as sfc
import scipy.fftpack as fft

# def analytic_m_calc(file,res):

scannum = 1760


path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
filespec = 'Msfixed'
cal_params_path = path+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

Ms = cal_params['Ms']
t = cal_params['t']
MstError = cal_params['MstError']
phi = cal_params['phi']
height = cal_params['height']
heightError = cal_params['heightError']

savepath = '/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/'

dwtypes = ["b"]
filespec = "test"

m = []
m.append(np.loadtxt(savepath+"mbx_"+filespec+".dat"))
m.append(np.loadtxt(savepath+"mby_"+filespec+".dat"))
m.append(np.loadtxt(savepath+"mbz_"+filespec+".dat"))

scd, vcd, meff, hk, h = sfc.stray_field_calc(m[0],m[1],m[2],Ms*t,2.5e-4,height)

np.savetxt(savepath+'b_x_'+filespec+'.txt', h[0], delimiter=',')
np.savetxt(savepath+'b_y_'+filespec+'.txt', h[1], delimiter=',')
np.savetxt(savepath+'b_z_'+filespec+'.txt', h[2], delimiter=',')

pi = np.pi
hzf = fft.fft2(h[2])
hzf = fft.fftshift(h[2])
scansize = 3.0e-4

dlen = len(hzf)
hlen = int(np.floor(dlen/2))
Vk = np.zeros_like(hzf)
k = 0
kmax = 2*pi*dlen/scansize

for j in range(0,dlen):
	ky = 2*pi*(j-hlen)/scansize
	for i in range(0,dlen):
		kx = 2*pi*(i-hlen)/scansize
		k = np.sqrt(kx**2 + ky**2)
		if (i==hlen and j==hlen):
			Vk[j,i] = 0
		else:
			Vk[j, i] = -hzf[j,i]/(k**2)

Vdata = np.real(fft.ifft2(fft.ifftshift(Vk)))

np.savetxt(savepath+'V_'+filespec+'.txt', Vdata, delimiter=',')
