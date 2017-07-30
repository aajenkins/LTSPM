# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T11:41:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-29T20:31:17-07:00



import numpy as np
import json
import matplotlib.pyplot as plt
import format_plot as fp
import stray_field_calc_fast as sfcf

def stray_field(scannum):

    path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
    scan_path = path+str(scannum)+'/'
    scan_params_path = scan_path+'scan_parameters.json'
    material_params_path = path+'material_parameters.json'
    with open(scan_params_path, 'r') as fread:
        scan_params = json.load(fread)
    with open(material_params_path, 'r') as fread:
        material_params = json.load(fread)

    Ms = material_params['Ms']
    MstError = material_params['MstError']
    t = material_params['t']

    phi = scan_params['phi']
    height = scan_params['height']
    heightError = scan_params['heightError']
    scanSize = scan_params['scanSize']
    zfield = scan_params['zfield']

    heights = [height - heightError, height, height + heightError]
    # MsList = [Ms - (MstError/t), Ms, Ms + (MstError/t)]
    Msts = [Ms*t - MstError, Ms*t, Ms*t + MstError]

    savepath = scan_path+'/stray_field_sim/'

    dwtypes = ["Bloch","RNeel", "LNeel"]
    errnames = ["lower", "mean", "upper"]

    for j in range(0,len(dwtypes)):
        filenames = ["mx"+dwtypes[j],"my"+dwtypes[j],"mz"]
        numfiles = len(filenames)
        m = []
        for i in range(0,numfiles):
            m.append(np.loadtxt(scan_path+filenames[i]+".txt", delimiter=','))

        for i in range(0,len(errnames)):
            print('calculating '+dwtypes[j]+' at '+errnames[i]+' height')

            h, scd, vcd, meff = sfcf.stray_field_calc_fast(m[0],m[1],m[2],Msts[i],scanSize,heights[i])

            h[2] = h[2] + zfield

            if (i==1):
                np.savetxt(savepath+'scd_'+dwtypes[j]+str(scannum)+'.txt', scd, delimiter=',')
                np.savetxt(savepath+'vcd_'+dwtypes[j]+str(scannum)+'.txt', vcd, delimiter=',')
                np.savetxt(savepath+'meff_'+dwtypes[j]+str(scannum)+'.txt', meff, delimiter=',')

            np.savetxt(savepath+dwtypes[j]+'_x_'+errnames[i]+'_'+str(scannum)+'.txt', h[0], delimiter=',')
            np.savetxt(savepath+dwtypes[j]+'_y_'+errnames[i]+'_'+str(scannum)+'.txt', h[1], delimiter=',')
            np.savetxt(savepath+dwtypes[j]+'_z_'+errnames[i]+'_'+str(scannum)+'.txt', h[2], delimiter=',')

if __name__ == "__main__":
    import sys
    if (len(sys.argv) == 2):
        stray_field(int(sys.argv[1]))
    else:
        print('enter scan number')
