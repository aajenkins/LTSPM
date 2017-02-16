# @Author: Jenkins Alec <alec>
# @Date:   2017-02-15T11:26:35-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-15T12:01:11-08:00

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
from scipy.ndimage.morphology import binary_erosion
from scipy.misc import imread

scannum = 1836
impath = '/Users/alec/UCSB/scan_images/irmn/domains'+str(scannum)+'contour.png'
domains = imread(impath, flatten=True)
bdomains = 1*domains
bdomains = np.multiply(1/255,domains).astype(int)

number_erosions = 15
image_sum = np.zeros((number_erosions))

for i in range(number_erosions):
    image_sum[i] = np.sum(bdomains)
    bdomains = binary_erosion(bdomains).astype(np.int)

plt.close('all')

fig1, ax1 = plt.subplots()
plt.imshow(bdomains)

fig1, ax1 = plt.subplots()
plt.plot(image_sum)

plt.show()
