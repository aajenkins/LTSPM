# @Author: Jenkins Alec <alec>
# @Date:   2017-02-10T13:40:16-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-12T15:17:35-08:00

import numpy as np
import m_reconstruction_Dovzhenko as mrecon
import time

t1 = time.time()

pi = np.pi

path = '/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/'
mzdata = np.loadtxt(path+'mbz_test.dat')

minmz = np.min(mzdata)
maxmz = np.max(mzdata)

mzdatanorm = np.multiply(np.add(mzdata,-(maxmz+minmz)/2),(2.0 - (1.0e-10))/(maxmz-minmz))
dlen = len(mzdatanorm);
phi_seed = np.full_like(mzdatanorm, 0)
for x in range(dlen):
    for y in range(dlen):
        phi_seed[x,y] = (pi/2)+np.arctan2(x-dlen/2,-(y-dlen/2))
np.savetxt(path+'mzdata_test.txt',mzdatanorm)

phi = mrecon.m_reconstruction_Dovzhenko(mzdatanorm, phi_seed, 50, 0.0001, 4*pi/180)

np.savetxt(path+'phi_seed_test.txt',phi_seed,delimiter=',')
np.savetxt(path+'phi_test.txt',phi,delimiter=',')

t2 = time.time()

print(t2-t1)
