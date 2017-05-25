# @Author: Jenkins Alec <alec>
# @Date:   2017-01-19T12:39:38-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-04-09T20:26:31-07:00



import numpy as np
import json
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc as sfc

# def analytic_m_calc(file,res):

scannum = 1903

path = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
filespec = 'Msfixed'
cal_params_path = path+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

Ms = cal_params['Ms']
t = cal_params['t']
MstError = cal_params['MstError']
phi = cal_params['phi']
scanSize = cal_params['scanSize']

height = 140.0e-9
dw0=15

basepath = "/Users/alec/UCSB/cofeb_analysis_data/irmn/stray_field_sim/"
savepath = "/Users/alec/UCSB/cofeb_analysis_data/irmn/stray_field_sim/heights/"

dwtypes = ['h', 'hr', 'nl', 'nr', 'b']
filenames = []

for i in range(0,len(dwtypes)):
    filenames.append(["m"+dwtypes[i]+"x_Msfixed",
                      "m"+dwtypes[i]+"y_Msfixed",
                      "m"+dwtypes[i]+"z_Msfixed"])

numdwtypes = len(dwtypes)
m = [[] for i in range(0,numdwtypes)]
for i in range(0,numdwtypes):
    for j in range(0,3):
        m[i].append(np.loadtxt(basepath+filenames[i][j]+".dat"))

for i in range(0, numdwtypes):
    print('calculating h for DW type '+str(dwtypes[i]))
    scd, vcd, meff, hk, h = sfc.stray_field_calc(m[i][0], m[i][1], m[i][2],
                                                 Ms*t, scanSize, height)

    slen = len(h[0])
    hlowres = [[], [], []]
    for k in range(0,3):
        hlowres[k] = h[k][0:slen:2, 0:slen:2]

    np.savetxt(savepath+dwtypes[i]+'_x_'+str(int(height*1e9))+'_'+str(scannum)+'.txt', hlowres[0], delimiter=',')
    np.savetxt(savepath+dwtypes[i]+'_y_'+str(int(height*1e9))+'_'+str(scannum)+'.txt', hlowres[1], delimiter=',')
    np.savetxt(savepath+dwtypes[i]+'_z_'+str(int(height*1e9))+'_'+str(scannum)+'.txt', hlowres[2], delimiter=',')
