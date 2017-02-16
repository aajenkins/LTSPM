# @Author: Jenkins Alec <alec>
# @Date:   2017-02-15T15:15:01-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-15T15:37:20-08:00

import scipy.fftpack as fft
import numpy as np


def calc_demag_energy_fourier(magn, res, thickness, Ms):

    slen = len(magn)
    hlen = slen/2

    surface_density_k = fft.fft2(magn)
    surface_density_k = fft.fftshift(surface_density_k)

    Hz_k = np.zeros_like(surface_density_k)
    Hz = np.zeros_like(surface_density_k)

    for j in range(0, slen):
        ky = 2*np.pi*(j-hlen)/(res*slen)
        for i in range(0,slen):
            kx = 2*np.pi*(i-hlen)/(res*slen)
            k = np.sqrt(kx**2 + ky**2)
            Hz_k[j][i] = k * surface_density_k[j][i]

    Hz_k = Ms * thickness * Hz_k
    Hz = np.real(fft.ifft2(fft.ifftshift(Hz_k)))

    energy = 0
    volume = res*res*thickness
    for j in range(0, slen):
        for i in range(0, slen):
            energy += Ms * volume * Hz[j][i] * magn[j][i]

    return energy, Hz
