# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 14:30:34 2017

@author: alec
"""

import numpy as np

def get_deltaC(mzsq, phisq, lphi):
	phiprime = phisq[1][1]+lphi
	dL=np.zeros((5))

	dL[0] = np.sqrt(1-(mzsq[1][1])**2)*((phisq[2][1]-phisq[0][1])*(np.cos(phiprime)-np.cos(phisq[1][1]))-
	(phisq[1][2]-phisq[1][0])*(np.sin(phiprime)-np.sin(phisq[1][1])))

	dL[1] = np.sqrt(1-(mzsq[1][2])**2)*lphi*np.sin(phisq[1][2])

	dL[2] = -1*np.sqrt(1-(mzsq[1][0])**2)*lphi*np.sin(phisq[1][0])

	dL[3] = -1*np.sqrt(1-(mzsq[1][0])**2)*lphi*np.cos(phisq[2][1])

	dL[4] = np.sqrt(1-(mzsq[1][0])**2)*lphi*np.cos(phisq[0][1])

	dLs = np.multiply(dL,dL)

	return np.sum(dLs)


def m_recon(mz, lphi, num_steps):

	phi = np.zeros_like(mz)
	deltaC = np.zeros_like(mz)
	l = len(mz)

	for k in range(0, num_steps):
		if (k%10 == 0):
			print('step '+str(k))
		for i in range(1,l-1):
			for j in range(1,l-1):
				deltaC[j][i] = get_deltaC(mz[j-1:j+2,i-1:i+2], phi[j-1:j+2,i-1:i+2], lphi)
				phi[j][i] = phi[j][i] - lphi*deltaC[j][i]

	return phi
