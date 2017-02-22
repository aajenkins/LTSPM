# @Author: Jenkins Alec <alec>
# @Date:   2017-02-17T13:16:12-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-17T13:19:12-08:00



import numpy as np
from scipy.ndimage import imread
import matplotlib.pyplot as plt
import format_plot as fp

scannum = 1836
impath = '/Users/alec/UCSB/scan_images/irmn/domains'+str(scannum)+'.png'
domains = imread(impath, flatten=True)
domains = np.multiply(1/255,domains)

dk = np.fft.fftshift(np.fft.fft2(domains))
drecon = np.real(np.fft.ifft2(np.fft.fftshift(dk)))


plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(domains, cmap='jet', interpolation='nearest')
plt.colorbar(im1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 0, 450)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(drecon, cmap='jet', interpolation='nearest')
plt.colorbar(im1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 550, 450)

plt.show()
