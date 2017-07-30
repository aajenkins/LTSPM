# @Author: Jenkins Alec <alec>
# @Date:   2017-02-28T14:33:58-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-27T16:30:15-07:00

import numpy as np

def calc_NV_field(bx, by, bz, theta, phi):
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+
    by*np.sin(theta)*np.sin(phi)+
    bz*np.cos(theta))

    return bnv

def calc_NV_field_angle(f1, f2, Dgs=2870, printNVCalcError=False):

    if (-(Dgs+f1+f2)*(Dgs+f1-2*f2)*(Dgs-2*f1+f2) < 0):
        Bnv = np.abs(f2-f1)/(2*2.8)
        Bp = 0
        if (printNVCalcError):
            print('Bnv imaginary',f1,f2)
    else:
        Bnv = (1/(3*np.sqrt(3*Dgs))) * np.sqrt( -(Dgs+f1+f2)*(Dgs+f1-2*f2)*(Dgs-2*f1+f2) )

        if (-(2*Dgs-f1-f2)*(2*Dgs-f1+2*f2)*(2*Dgs+2*f1-f2) < 0):
            Bnv = np.abs(f2-f1)/(2*2.8)
            Bp = 0
            if (printNVCalcError):
                print('Bp imaginary',f1,f2)
        else:
            Bp = (1/(3*np.sqrt(3*Dgs))) * np.sqrt( -(2*Dgs-f1-f2)*(2*Dgs-f1+2*f2)*(2*Dgs+2*f1-f2) )

    theta = np.arctan2(Bp, Bnv)

    return Bnv, Bp, theta
