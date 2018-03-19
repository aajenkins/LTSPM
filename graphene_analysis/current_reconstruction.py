# @Author: Jenkins Alec <alec>
# @Date:   2018-03-15T22:09:25-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2018-03-16T21:59:14-07:00

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np

#---------------- RECONSTRUCTION ---------------------------------
#-----------------------------------------------------------------

def current_reconstruction_1D(bnv, theta, phi, z, scansize):

    pi = np.pi
    bnvf = np.fft.fft(bnv)
    bnvf = np.fft.fftshift(bnvf)

    dlen = len(bnvf)
    hlen = int(np.floor(dlen/2))

    jyf = np.zeros_like(bnvf)

    for i in range(0,dlen):
        kx = 2*pi*(i-hlen)/scansize
        if (i == hlen):
            jyf = 0
        else:
            jyf = bnvf/( (2e-7) * np.exp(-z*np.abs(kx)) * (np.sin(theta)*np.cos(phi)
                                                           - 1j*np.sign(kx)*np.cos(theta)) )

    jy = np.real(np.fft.ifft(np.fft.ifftshift(jyf)))

    return jy, jyf

def current_reconstruction_afm_1D(bnv, theta, phi, afm, scansize):

    pi = np.pi
    bnvf = np.fft.fft(bnv)
    bnvf = np.fft.fftshift(bnvf)

    dlen = len(bnvf)
    hlen = int(np.floor(dlen/2))

    jyf = np.zeros_like(bnvf)

    for i in range(0,dlen):
        kx = 2*pi*(i-hlen)/scansize
        if (i == hlen):
            jyf = 0
        else:
            jyf = bnvf/( (2e-7) * np.exp(-z*np.abs(kx)) * (np.sin(theta)*np.cos(phi)
                                                           - 1j*np.sign(kx)*np.cos(theta)) )

    jy = np.real(np.fft.ifft(np.fft.ifftshift(jyf)))

    return jy, jyf
