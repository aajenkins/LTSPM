# @Author: Jenkins Alec <alec>
# @Date:   2017-02-03T18:34:58-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-04-07T10:51:20-07:00

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import imread
import json

impath = '/Users/alec/UCSB/scan_images/irmn/domains1836.png'
domains = imread(impath, flatten=True)

path = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
filespec = 'Msfixed'
cal_params_path = path+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

scansize = 2*(5e-6)
slen = len(domains)
res = scansize/slen
res_difference = 2e-9
thickness = cal_params['t']
Ms = cal_params['Ms']

domains = np.add(np.multiply(2/255,domains),-1)

ovfPath = '/Users/alec/UCSB/cofeb_analysis_data/irmn/oommf/'
ovfHeaderFile = open(ovfPath+'irmn_magn_header.ovf', 'r')
ovfHeaderFileDres = open(ovfPath+'irmn_magn_header_dres.ovf', 'r')

ovfDomains = np.zeros((slen**2,3))

for j in range(slen):
    for i in range(slen):
        ovfDomains[slen*j+i, 2] = domains[j][i]

ovfHeader = ''
ovfHeaderDres = ''
for l in ovfHeaderFile.readlines():
    ovfHeader += str(l)
for l in ovfHeaderFileDres.readlines():
    ovfHeaderDres += str(l)
# ovfHeader += str(ovfHeaderFile.readlines()[:-1])
# ovfHeaderDres += str(ovfHeaderFileDres.readlines()[:-1])

np.savetxt(ovfPath+'irmn_magn.ovf', ovfDomains, header=ovfHeader, comments='')
np.savetxt(ovfPath+'irmn_magn_dres.ovf', ovfDomains, header=ovfHeaderDres, comments='')
ovfHeaderFile.close()
ovfHeaderFileDres.close()

ovfFile = open(ovfPath+'irmn_magn.ovf', 'a')
ovfFile.write('# End: data text \n# End: segment')
ovfFile.close()

ovfFileDres = open(ovfPath+'irmn_magn_dres.ovf', 'a')
ovfFileDres.write('# End: data text \n# End: segment')
ovfFileDres.close()
