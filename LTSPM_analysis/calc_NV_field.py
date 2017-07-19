# @Author: Jenkins Alec <alec>
# @Date:   2017-02-28T14:33:58-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-19T09:27:57-07:00

import numpy as np

def calc_NV_field(bx, by, bz, theta, phi):
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+
    by*np.sin(theta)*np.sin(phi)+
    bz*np.cos(theta))

    return bnv

def calc_NV_field_angle(f1, f2, Dgs=2870):
    if ((f1+f2)/2 < Dgs-(1e-10)):
        B0 = 1000
        Bnv = B0
        theta = 0
        print('Bnv< fail',f1,f2)
    # elif (Dgs-1 < (f1+f2)/2 < Dgs):
    #     Bnv = np.abs(f1-f2)/(2*2.8)
    #     B0 = Bnv
    #     theta = 0
    else:
        P = f1**2 + f2**2 - f1*f2
        Q = (f1+f2)*(2*(f1**2) + 2*(f2**2) - 5*f1*f2)
        if ((P - (Dgs**2)) < 0):
            B0 = 1000
            Bnv = B0
            theta = 0
            print('P< fail',f1,f2)
        else:
            B0 = np.sqrt((P - (Dgs**2))/3)
            cos_sq = (Q + 9*Dgs*(B0**2) + 2*(Dgs**3)) / (27*Dgs*(B0**2))
            if (cos_sq > 1):
                print('cos fail', cos_sq)
                cos_sq=1
            if (cos_sq < 0):
                print('cos fail', cos_sq)
                cos_sq=0
            theta = np.arccos(np.sqrt( cos_sq ))
            Bnv = B0*np.cos(theta)

    return Bnv, B0, theta
