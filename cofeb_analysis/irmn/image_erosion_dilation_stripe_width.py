# @Author: Jenkins Alec <alec>
# @Date:   2017-02-15T11:26:35-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-15T14:42:57-08:00

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
from scipy.ndimage.morphology import binary_erosion, binary_dilation
from scipy.misc import imread
from scipy.optimize import curve_fit

def rolloff(x, *args):
    imsum = np.zeros_like(x)
    c = args[0]
    x0 = args[1]
    n=2
    imsum = c/(x**n + x0**n)
    return imsum

scannum = 1836
impath = '/Users/alec/UCSB/scan_images/irmn/domains'+str(scannum)+'contour.png'
domains = imread(impath, flatten=True)
erosion_domains = 1*domains
erosion_domains = np.multiply(1/255,domains).astype(int)
dilation_domains = 1*domains
dilation_domains = np.multiply(1/255,domains).astype(int)

number_steps = 15
image_sum_erosion = np.zeros((number_steps))
image_sum_dilation = np.zeros((number_steps))

for i in range(number_steps):
    image_sum_erosion[i] = np.sum(erosion_domains)
    erosion_domains = binary_erosion(erosion_domains).astype(np.int)
    image_sum_dilation[i] = np.sum(dilation_domains)
    dilation_domains = binary_dilation(dilation_domains).astype(np.int)
image_sum_dilation = image_sum_dilation[-1]-image_sum_dilation
#
# x = np.arange(number_steps)
# xfine = np.arange(0,number_steps,0.2)
# guess = [1.0e4, 5.5]
# popt, pcov = curve_fit(rolloff, x, image_sum_erosion, p0=guess)

erosion_slope = image_sum_erosion[1]-image_sum_erosion[0]
erosion_line = erosion_slope * np.arange(8)+image_sum_erosion[0]
erosion_intercept = -image_sum_erosion[0]/erosion_slope

dilation_slope = image_sum_dilation[1]-image_sum_dilation[0]
dilation_line = dilation_slope * np.arange(8)+image_sum_dilation[0]
dilation_intercept = -image_sum_dilation[0]/dilation_slope

mean_pixel_width = 2*((erosion_intercept + dilation_intercept)/2)

print('erosion x axis intercept = '+str(erosion_intercept))
print('dilation x axis intercept = '+str(dilation_intercept))

print('mean +/- pixel width = '+str(mean_pixel_width))

plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(image_sum_dilation)
plt.plot(dilation_line,'r-')
plt.ylim([0,1.2e4])

fig1, ax1 = plt.subplots()
plt.plot(image_sum_erosion)
plt.plot(erosion_line,'r-')
plt.ylim([0,1.2e4])

plt.show()
