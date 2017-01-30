# @Author: Jenkins Alec <alec>
# @Date:   2017-01-25T16:16:05-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-01-29T21:10:56-08:00

# calculate the demagnetization energy from a stripe phase contour scan

import numpy as np
from scipy import signal

pi = np.pi
mu0 = 4*pi*(1e-7)

def f(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    f = (r/6)*(3*(z**2) - x**2 - y**2)
    if (y > 0):
        f += 0.5 * x * (y**2 - z**2)*np.arcsinh(x/(y**2 + z**2))
    if (x > 0):
        f += 0.5 * y * (x**2 - z**2)*np.arcsinh(y/(x**2 + z**2))
    if (z>0):
        f -= x*y*z*np.arctan(x*y/(z*r))
    return f

def F2(x,y,z):
    return ( f(x,y,z) - f(0,y,z) - f(x,0,z) + f(0,0,z) )

def F1(x,y,z,a):
    return ( F2(x,y,z) - F2(x-a,y,z) - F2(x,y-a,z) - F2(x-a,y-a,z) )

def F(x,y,z,a):
    return ( F1(x+a,y+a,z,a) - F1(x,y+a,z,a) - F1(x+a,y,z,a) + F1(x,y,z,a) )

def Nzz(x,y,z,a,t):
    return (1/(4*pi))*( 2*F(x,y,z,a) - F(x,y,z+t,a) - F(x,y,z-t,a) )

def get_demag_energy(magn, Ms, a, t):
    slen = len(magn)
    demag_tensor = np.zeros_like(magn)

    for j in range(0, slen):
        for i in range(0, slen):
            demag_tensor[j,i] = Nzz(i*a, j*a, 0, a, t)

    h = signal.convolve2d(demag_tensor, magn, boundary='symm')
    energy = 0

    for j in range(0, slen):
        for i in range(0, slen):
            energy += h[j,i]*magn[j,i]

    np.multiply(energy, -(a**2)*(Ms**2)*mu0/2)

    return energy
