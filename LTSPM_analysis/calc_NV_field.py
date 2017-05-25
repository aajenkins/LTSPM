# @Author: Jenkins Alec <alec>
# @Date:   2017-02-28T14:33:58-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-28T15:40:44-08:00

import numpy as np

def calc_NV_field(bx, by, bz, theta, phi):
    bnv = np.abs(bx*np.sin(theta)*np.cos(phi)+
    by*np.sin(theta)*np.sin(phi)+
    bz*np.cos(theta))

    return bnv
