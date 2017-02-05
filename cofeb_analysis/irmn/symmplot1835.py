# @Author: Jenkins Alec <alec>
# @Date:   2017-01-31T21:02:36-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-03T15:11:35-08:00



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
domains1 = np.flipud(np.fliplr(domains))
domains2 = np.flipud(domains)
domains3 = np.flipud(np.fliplr(domains))
domains4 = np.fliplr(domains)
domains5 = np.fliplr(domains)
domains6 = np.flipud(np.fliplr(domains))
domains7 = np.flipud(domains)
domains8 = np.flipud(np.fliplr(domains))

domainstop = np.concatenate((domains1,domains2,domains3), axis=1)
domainsmiddle = np.concatenate((domains4,domains,domains5), axis=1)
domainsbottom = np.concatenate((domains6,domains7,domains8), axis=1)

domains_mirror = np.concatenate((domainstop, domainsmiddle, domainsbottom),
                                axis=0)



plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(domains_mirror, cmap='Greys', interpolation='nearest')
fp.format_plot(plt, 800, 800, 50, 50)

plt.show()
