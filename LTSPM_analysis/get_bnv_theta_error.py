# @Author: Jenkins Alec <alec>
# @Date:   2017-02-28T15:24:12-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-28T15:28:07-08:00



import numpy as np

def get_bnv_theta_error(bx,by,bz,theta,thetaError,phi):
    bnvError = np.abs( bx*np.cos(theta)*np.cos(phi)
                      + by*np.cos(theta)*np.sin(phi)
                      - bz*np.sin(theta) )

    return bnvError
