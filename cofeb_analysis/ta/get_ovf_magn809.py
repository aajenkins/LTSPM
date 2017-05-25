# @Author: Jenkins Alec <alec>
# @Date:   2017-02-03T18:34:58-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-04-10T11:09:51-07:00

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import imread

impath = '/Users/alec/UCSB/scan_images/ta_contour_selection/domains809.png'
domains = imread(impath, flatten=True)

slen = len(domains)

domains = np.add(np.multiply(2/255,domains),-1)

ovfPath = '/Users/alec/UCSB/oommf/data_and_runs/ta/'
ovfHeaderFile1 = open(ovfPath+'ta_magn_header1.ovf', 'r')
ovfHeaderFileDres1 = open(ovfPath+'ta_magn_header_dres1.ovf', 'r')
ovfHeaderFile2 = open(ovfPath+'ta_magn_header2.ovf', 'r')
ovfHeaderFileDres2 = open(ovfPath+'ta_magn_header_dres2.ovf', 'r')

ovfDomains = np.zeros((slen**2,3))

for j in range(slen):
    for i in range(slen):
        ovfDomains[slen*j+i, 2] = domains[j][i]

ovfHeader1 = ''
ovfHeaderDres1 = ''
ovfHeader2 = ''
ovfHeaderDres2 = ''
for l in ovfHeaderFile1.readlines():
    ovfHeader1 += str(l)
for l in ovfHeaderFileDres1.readlines():
    ovfHeaderDres1 += str(l)
for l in ovfHeaderFile2.readlines():
    ovfHeader2 += str(l)
for l in ovfHeaderFileDres2.readlines():
    ovfHeaderDres2 += str(l)
# ovfHeader += str(ovfHeaderFile.readlines()[:-1])
# ovfHeaderDres += str(ovfHeaderFileDres.readlines()[:-1])

np.savetxt(ovfPath+'ta_magn1.ovf', ovfDomains, header=ovfHeader1, comments='')
np.savetxt(ovfPath+'ta_magn_dres1.ovf', ovfDomains, header=ovfHeaderDres1, comments='')
np.savetxt(ovfPath+'ta_magn2.ovf', ovfDomains, header=ovfHeader2, comments='')
np.savetxt(ovfPath+'ta_magn_dres2.ovf', ovfDomains, header=ovfHeaderDres2, comments='')
ovfHeaderFile1.close()
ovfHeaderFileDres1.close()
ovfHeaderFile2.close()
ovfHeaderFileDres2.close()

ovfFile1 = open(ovfPath+'ta_magn1.ovf', 'a')
ovfFile1.write('# End: data text \n# End: segment')
ovfFile1.close()

ovfFileDres1 = open(ovfPath+'ta_magn_dres1.ovf', 'a')
ovfFileDres1.write('# End: data text \n# End: segment')
ovfFileDres1.close()

ovfFile2 = open(ovfPath+'ta_magn2.ovf', 'a')
ovfFile2.write('# End: data text \n# End: segment')
ovfFile2.close()

ovfFileDres2 = open(ovfPath+'ta_magn_dres2.ovf', 'a')
ovfFileDres2.write('# End: data text \n# End: segment')
ovfFileDres2.close()
