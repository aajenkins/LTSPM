# @Author: Jenkins Alec <alec>
# @Date:   2017-02-15T15:15:01-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-17T13:38:05-08:00

import numpy as np


def calc_demag_energy_fourier(magn, res, thickness, Ms):

    slen = len(magn)
    hlen = slen/2
    Mst = Ms*thickness
    sprefactor = (1e-3)*Mst

    surface_density_k = np.fft.fftshift(np.fft.fft2(magn, norm="ortho"))

    Hz_k = np.zeros_like(surface_density_k)
    Hz = np.zeros_like(surface_density_k)

    for j in range(0, slen):
        ky = 2*np.pi*(j-hlen)/(res*slen)
        for i in range(0,slen):
            kx = 2*np.pi*(i-hlen)/(res*slen)
            k = np.sqrt(kx**2 + ky**2)
            Hz_k[j][i] = 2*np.pi * k * surface_density_k[j][i]

    Hz = sprefactor * np.real(np.fft.ifft2(np.fft.ifftshift(Hz_k)))
    # Hzb = Hz - (np.max(Hz)-np.min(Hz))/2

    energy = 0
    volume = res*res*thickness
    for j in range(0, slen):
        for i in range(0, slen):
            energy += Hz[j][i] * magn[j][i]

    energy = Ms * volume * energy
    return energy, Hz, Hz_k
