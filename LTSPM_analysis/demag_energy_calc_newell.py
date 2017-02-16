# @Author: Jenkins Alec <alec>
# @Date:   2017-01-25T16:16:05-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-14T20:43:08-08:00

# calculate the demagnetization energy from a stripe phase contour scan

import numpy as np
from scipy import signal

import symmetric_padding_image as spi

pi = np.pi
mu0 = 4*pi*(1e-7)

def demag_energy_calc_newell(magn, demag_tensor, Ms, a, t):

    slen = len(magn)
    magn_mirror = spi.symmetric_padding_image(magn)
    magn_mirror = np.multiply(magn_mirror, Ms)

    h = Ms * signal.convolve2d(demag_tensor, magn_mirror, mode='valid')
    energy = 0

    h = h[2:slen,2:slen]

    for j in range(0, slen-2):
        for i in range(0, slen-2):
            energy += h[j,i]*magn[j+1,i+1]

    energy = -(a*a*t)*(Ms)*(mu0/2)*energy

    return energy, h
