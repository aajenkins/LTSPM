# @Author: Jenkins Alec <alec>
# @Date:   2017-02-12T13:28:44-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-12T13:37:57-08:00

import numpy as np

pi = np.pi

phi1 = -(2*pi-0.01)
phi0 = 0

if ((phi1-phi0)%(2*pi) > pi):
    d_phi = (phi1-phi0)%(2*pi)-2*pi
else:
    d_phi = (phi1-phi0)%(2*pi)

print(d_phi)
