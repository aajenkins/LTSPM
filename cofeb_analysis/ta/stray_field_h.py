# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T11:41:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-10T13:16:24-07:00



import numpy as np
import json
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc as sfc

def stray_field_h(scannum):

    path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
    cal_params_path = path+'cal_parameters.json'
    with open(cal_params_path, 'r') as fread:
        cal_params = json.load(fread)

    Ms = cal_params['Ms']
    t = cal_params['t']
    MstError = cal_params['MstError']
    phi = cal_params['phi']
    height = cal_params['height']
    heightError = cal_params['heightError']
    scanSize = cal_params['scanSize']

    heights = [height - heightError, height, height + heightError]
    Msts = [Ms*t - MstError, Ms*t, Ms*t + MstError]

    savepath = '/Users/alec/UCSB/cofeb_analysis_data/ta/'+str(scannum)+'/stray_field_sim/'

    dwtypes = ["Bloch","RNeel", "LNeel"]
    errnames = ["lower", "mean", "upper"]

    for j in range(0,len(dwtypes)):
        filenames = ["mx"+dwtypes[j],"my"+dwtypes[j],"mz"]
        numfiles = len(filenames)
        m = []
        for i in range(0,numfiles):
            m.append(np.loadtxt(savepath+filenames[i]+".dat"))

        for i in range(0,len(errnames)):
            print('calculating '+dwtypes[j]+' at '+errnames[i]+' height')

            scd, vcd, meff, hk, h = sfc.stray_field_calc(m[0],m[1],m[2],Msts[i],scanSize,heights[i])

            np.savetxt(savepath+dwtypes[j]+'_x_'+errnames[i]+'_'+str(scannum)+'.txt', h[0], delimiter=',')
            np.savetxt(savepath+dwtypes[j]+'_y_'+errnames[i]+'_'+str(scannum)+'.txt', h[1], delimiter=',')
            np.savetxt(savepath+dwtypes[j]+'_z_'+errnames[i]+'_'+str(scannum)+'.txt', h[2], delimiter=',')
