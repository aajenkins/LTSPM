# @Author: Jenkins Alec <alec>
# @Date:   2017-02-28T14:33:58-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-19T20:13:57-07:00

import numpy as np

def calc_NV_field(bx, by, bz, theta, phi):
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+
    by*np.sin(theta)*np.sin(phi)+
    bz*np.cos(theta))

    return bnv

def calc_NV_field_angle(f1, f2, Dgs=2870):

    if (-(Dgs+f1+f2)*(Dgs+f1-2*f2)*(Dgs-2*f1+f2) < 0):
        Bnv = np.abs(f2-f1)/(2*2.8)
        Bp = 0
        print('Bnv imaginary',f1,f2)
    # elif (Dgs-1 < (f1+f2)/2 < Dgs):
    #     Bnv = np.abs(f1-f2)/(2*2.8)
    #     B0 = Bnv
    #     theta = 0
    else:
        Bnv = (1/(3*np.sqrt(3*Dgs))) * np.sqrt( -(Dgs+f1+f2)*(Dgs+f1-2*f2)*(Dgs-2*f1+f2) )

        if (-(2*Dgs-f1-f2)*(2*Dgs-f1+2*f2)*(2*Dgs+2*f1-f2) < 0):
            Bnv = np.abs(f2-f1)/(2*2.8)
            Bp = 0
            print('Bp imaginary',f1,f2)
        else:
            Bp = (1/(3*np.sqrt(3*Dgs))) * np.sqrt( -(2*Dgs-f1-f2)*(2*Dgs-f1+2*f2)*(2*Dgs+2*f1-f2) )

    tantheta = Bp/Bnv
    theta = np.arctan(tantheta)

    return Bnv, Bp, theta
