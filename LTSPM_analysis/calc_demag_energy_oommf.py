# @Author: Jenkins Alec <alec>
# @Date:   2017-02-17T14:20:47-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-17T14:23:22-08:00



import numpy as np

ohfPath = '/Users/alec/UCSB/oommf/data_and_runs/irmn/'
filename = 'irmn_magn-hdemag.ohf'
ohfDemagH = np.genfromtxt(ohfPath+filename, skip_header=32, skip_footer=2)

print(ohfDemagH[0])
