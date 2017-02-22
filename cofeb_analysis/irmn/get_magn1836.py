# @Author: Jenkins Alec <alec>
# @Date:   2017-02-03T18:34:58-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-17T15:26:41-08:00

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import imread

impath = '/Users/alec/UCSB/scan_images/irmn/domains1836.png'
domains = imread(impath, flatten=True)

scanvsize = 2
scansize = (5e-6)*scanvsize
slen = len(domains)
res = scansize/slen
res_difference = 2e-9
thickness = 0.911e-9
Ms = 1.178e6

domains = np.add(np.multiply(2/255,domains),-1)

ovfPath = '/Users/alec/UCSB/cofeb_analysis_data/irmn/oommf/'
ovfheaderfile = open(ovfPath+'irmn_magn_header_dres.ovf', 'r')

domainsovf = np.zeros((slen**2,3))

for j in range(slen):
    for i in range(slen):
        domainsovf[slen*j+i, 2] = domains[j][i]

np.savetxt(ovfPath+'irmn_magn_dres.ovf', domainsovf, header=ovfheaderfile.read(), comments='')
ovfheaderfile.close()

ovffile = open(ovfPath+'irmn_magn_dres.ovf', 'a')
ovffile.write('# End: data text \n# End: segment')
ovffile.close()
