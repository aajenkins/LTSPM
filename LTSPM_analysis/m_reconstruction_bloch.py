# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T12:05:35-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-16T08:03:36-07:00

import numpy as np

pi = np.pi


def m_reconstruction(mz):
	mzlen = len(mz)

	mzgrad = np.gradient(mz)

	mzgradx = mzgrad[1]
	mzgrady = mzgrad[0]

	phi_bloch = pi/2 + np.arctan2(mzgrady, mzgradx)

	mx = np.sqrt(1-mz**2)*np.cos(phi_bloch)
	my = np.sqrt(1-mz**2)*np.sin(phi_bloch)

	return phi_bloch, mx, my
