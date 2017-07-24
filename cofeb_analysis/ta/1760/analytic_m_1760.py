# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T11:41:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-22T14:34:57-07:00



import numpy as np
import json
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc_fast as sfcf

# def analytic_m_calc(file,res):

scannum = 1760


path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
scan_params_path = path+str(scannum)+'/'+'scan_parameters.json'
material_params_path = path+'material_parameters.json'
with open(scan_params_path, 'r') as fread:
    scan_params = json.load(fread)
with open(material_params_path, 'r') as fread:
    material_params = json.load(fread)

Ms = material_params['Ms']
t = material_params['t']
MstError = material_params['MstError']
phi = scan_params['phi']
height = scan_params['height']
heightError = scan_params['heightError']
simSize = scan_params['scanSize']

heights = [height - heightError, height, height + heightError]
Msts = [Ms*t - MstError, Ms*t, Ms*t + MstError]

savepath = '/Users/alec/UCSB/cofeb_analysis_data/ta/1760/stray_field_sim/'

# dwtypes = ["Bloch", "LNeel", "RNeel"]
dwtypes = ["Bloch"]
errnames = ["lower", "mean", "upper"]

for j in range(0,len(dwtypes)):
    filenames = ["mx"+dwtypes[j],"my"+dwtypes[j],"mz"]
    numfiles = len(filenames)
    m = []
    for i in range(0,numfiles):
        m.append(np.loadtxt(savepath+filenames[i]+".txt", delimiter=','))

    for i in range(0,len(errnames)):
        print('calculating '+dwtypes[j]+' at '+errnames[i]+' height')

        scd, vcd, meff, hk, h = sfcf.stray_field_calc_fast(m[0],m[1],m[2],Msts[i],simSize,heights[i])

        # np.savetxt(savepath+dwtypes[j]+'_x_'+errnames[i]+'_'+str(scannum)+filespec+'.txt', h[0], delimiter=',')
        # np.savetxt(savepath+dwtypes[j]+'_y_'+errnames[i]+'_'+str(scannum)+filespec+'.txt', h[1], delimiter=',')
        # np.savetxt(savepath+dwtypes[j]+'_z_'+errnames[i]+'_'+str(scannum)+filespec+'.txt', h[2], delimiter=',')

        slen = len(h[0])
        hlowres = [[], [], []]
        for k in range(0,numfiles):
            hlowres[k] = h[k][0:slen:2, 0:slen:2]

        np.savetxt(savepath+dwtypes[j]+'_x_'+errnames[i]+'_lowres.txt', hlowres[0], delimiter=',')
        np.savetxt(savepath+dwtypes[j]+'_y_'+errnames[i]+'_lowres.txt', hlowres[1], delimiter=',')
        np.savetxt(savepath+dwtypes[j]+'_z_'+errnames[i]+'_lowres.txt', hlowres[2], delimiter=',')
