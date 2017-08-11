# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T11:41:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-10T09:50:44-07:00



import numpy as np
import json
import stray_field_calc_thick as sfct

pi = np.pi


def stray_field(scannum, helicities = [0, 90, 180], errnames = ["lower", "mean", "upper"]):

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
    Mss = [Ms - MstError/t, Ms, Ms + MstError/t]

    savepath = scan_path+'stray_field_sim/'

    for j in range(0,len(helicities)):
        filenames = ["mx"+str(helicities[j]),"my"+str(helicities[j]),"mz"]
        numfiles = len(filenames)
        m = []
        for i in range(0,numfiles):
            m.append(np.loadtxt(scan_path+filenames[i]+".txt", delimiter=','))

        for i in range(0,len(errnames)):
            print('calculating B for helicity '+str(helicities[j])+' at '+errnames[i]+' height')

            h, scd, vcd, meff = sfct.stray_field_calc_thick(m[0], m[1], m[2], Mss[i], t, scanSize, heights[i])

            h[2] = h[2] + zfield

            slen = len(h[0])
            if(slen > 200):
                scaleFactor = int(slen/200)
                h[0] = h[0][::scaleFactor, ::scaleFactor]
                h[1] = h[1][::scaleFactor, ::scaleFactor]
                h[2] = h[2][::scaleFactor, ::scaleFactor]
                scd = scd[::scaleFactor, ::scaleFactor]
                vcd = vcd[::scaleFactor, ::scaleFactor]
                meff = meff[::scaleFactor, ::scaleFactor]

            if (i==1):
                np.savetxt(savepath+'scd_'+str(helicities[j])+str(scannum)+'.txt', scd, delimiter=',')
                np.savetxt(savepath+'vcd_'+str(helicities[j])+str(scannum)+'.txt', vcd, delimiter=',')
                np.savetxt(savepath+'meff_'+str(helicities[j])+str(scannum)+'.txt', meff, delimiter=',')

            np.savetxt(savepath+'h'+str(helicities[j])+'_x_'+errnames[i]+'_'+str(scannum)+'.txt', h[0], delimiter=',')
            np.savetxt(savepath+'h'+str(helicities[j])+'_y_'+errnames[i]+'_'+str(scannum)+'.txt', h[1], delimiter=',')
            np.savetxt(savepath+'h'+str(helicities[j])+'_z_'+errnames[i]+'_'+str(scannum)+'.txt', h[2], delimiter=',')

if __name__ == "__main__":
    import sys
    if (len(sys.argv) == 2):
        stray_field(int(sys.argv[1]))
    elif (len(sys.argv) == 3):
        stray_field(int(sys.argv[1]), helicities=np.array(eval(sys.argv[2])))
    else:
        print('enter scan number')
