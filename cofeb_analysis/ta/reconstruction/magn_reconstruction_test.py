# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:56:50 2016

@author: alec
"""

#---------------- IMPORTS ----------------------------------------
#-----------------------------------------------------------------

import numpy as np
import m_reconstruction_Dovzhenko as mr

pi = np.pi

path = '/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/'
mzdata = np.loadtxt(path+'mbz_test.dat')

minmz = np.min(mzdata)
maxmz = np.max(mzdata)

mzdatanorm = np.multiply(np.add(mzdata,-(maxmz+minmz)/2),(2.0 - (1.0e-10))/(maxmz-minmz))

np.savetxt(path+'mzdata_test.txt',mzdatanorm,delimiter=',')

phi = mr.m_reconstruction_Dovzhenko(mzdatanorm, 20, 0.001, 1*pi/180)

np.savetxt(path+'phi_test.txt',phi,delimiter=',')
