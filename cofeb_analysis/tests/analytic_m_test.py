# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T11:41:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-03-17T21:01:19-05:00



import numpy as np
import json
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc as sfc

# def analytic_m_calc(file,res):

scannum = 1760

calpath = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
filespec = 'Msfixed'
cal_params_path = calpath+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

Ms = cal_params['Ms']
t = cal_params['t']
MstError = cal_params['MstError']
phi = cal_params['phi']
height = cal_params['height']
heightError = cal_params['heightError']

heights = [height - heightError, height, height + heightError]
Msts = [Ms*t - MstError, Ms*t, Ms*t + MstError]

savepath = '/Users/alec/UCSB/mathematica/CoFeB-MgO/skyrmion_shape/tests/'

filenames = ['mnrx_test.dat', 'mnry_test.dat', 'mnrz_test.dat']

numfiles = len(filenames)
m = []
for i in range(0,numfiles):
    m.append(np.loadtxt(savepath+filenames[i]))

scd, vcd, meff, hk, h = sfc.stray_field_calc(m[0],m[1],m[2],Ms*t,2.5e-4,height)

slen = len(h[0])
hlowres = [[], [], []]
for k in range(0,numfiles):
    hlowres[k] = h[k][0:slen:10, 0:slen:10]

np.savetxt(savepath+'b_x.txt', hlowres[0], delimiter=',')
np.savetxt(savepath+'b_y.txt', hlowres[1], delimiter=',')
np.savetxt(savepath+'b_z.txt', hlowres[2], delimiter=',')
