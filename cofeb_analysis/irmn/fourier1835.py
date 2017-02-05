# @Author: Jenkins Alec <alec>
# @Date:   2017-01-31T21:02:36-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-01-31T21:14:21-08:00



import matplotlib.pyplot as plt
import numpy as np
import scipy.fftpack as fft
from scipy.ndimage import imread

import format_plot as fp
import fourier_image as fi

impath = '/Users/alec/UCSB/scan_images/irmn/domains1835.png'
domains = imread(impath, flatten=True)

scanvsize = 2
scansize = (5e-6)*scanvsize
slen = len(domains)
res = scansize/slen
res_difference = 2e-9
thickness = 0.911e-9
Ms = 1.044e6


domains = np.add(np.multiply(2/255,domains),-1)

wdomains = fi.window_image(domains)
fdomains = fft.fft2(wdomains)
fdomains = fft.fftshift(fdomains)


plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(domains, cmap='Greys', interpolation='nearest')
fp.format_plot(plt, 400, 400, 450, 50)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(np.log(np.abs(fdomains)+1), cmap='Greys', interpolation='nearest')
fp.format_plot(plt, 400, 400, 450, 50)

plt.show()
