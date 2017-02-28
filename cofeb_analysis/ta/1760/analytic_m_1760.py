# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T11:41:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-28T14:18:17-08:00



import numpy as np
import json
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc as sfc

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


heights = [height - heightError, height, height + heightError]
Msts = [Ms*t - MstError, Ms*t, Ms*t + MstError]

savepath = '/Users/alec/UCSB/cofeb_analysis_data/ta/stray_field_sim/'

dwtypes = ["nr", "nl", "b"]
errnames = ["lower", "mean", "upper"]
filespec = "Msfixed"

for j in range(0,len(dwtypes)):
    filenames = ["m"+dwtypes[j]+"x","m"+dwtypes[j]+"y","m"+dwtypes[j]+"z"]
    numfiles = len(filenames)
    m = []
    for i in range(0,numfiles):
        m.append(np.loadtxt(savepath+filenames[i]+"_"+filespec+".dat"))

    for i in range(0,len(errnames)):
        print('calculating '+dwtypes[j]+' at '+errnames[i]+' height')

        scd, vcd, meff, hk, h = sfc.stray_field_calc(m[0],m[1],m[2],Msts[i],2.5e-4,heights[i])

        np.savetxt(savepath+dwtypes[j]+'_x_'+errnames[i]+'_'+str(scannum)+filespec+'.txt', h[0], delimiter=',')
        np.savetxt(savepath+dwtypes[j]+'_y_'+errnames[i]+'_'+str(scannum)+filespec+'.txt', h[1], delimiter=',')
        np.savetxt(savepath+dwtypes[j]+'_z_'+errnames[i]+'_'+str(scannum)+filespec+'.txt', h[2], delimiter=',')

        slen = len(h[0])
        hlowres = [[], [], []]
        for k in range(0,numfiles):
            hlowres[k] = h[k][0:slen:2, 0:slen:2]

        np.savetxt(savepath+dwtypes[j]+'_x_'+errnames[i]+'_lowres_'+str(scannum)+filespec+'.txt', hlowres[0], delimiter=',')
        np.savetxt(savepath+dwtypes[j]+'_y_'+errnames[i]+'_lowres_'+str(scannum)+filespec+'.txt', hlowres[1], delimiter=',')
        np.savetxt(savepath+dwtypes[j]+'_z_'+errnames[i]+'_lowres_'+str(scannum)+filespec+'.txt', hlowres[2], delimiter=',')
