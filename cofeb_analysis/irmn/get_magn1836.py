# @Author: Jenkins Alec <alec>
# @Date:   2017-02-03T18:34:58-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-04T20:54:57-08:00



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
Ms = 1.044e6

domains = np.add(np.multiply(2/255,domains),-1)
