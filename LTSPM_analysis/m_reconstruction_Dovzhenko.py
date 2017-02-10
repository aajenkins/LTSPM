# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T12:05:35-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-09T17:01:04-08:00



import numpy as np

def lagr(mz, phi, x, y):
	mr = np.sqrt(1 - (mz[x,y])**2)
	dxphi = (phi[x+1,y]-phi[x-1,y])/2
	dyphi = (phi[x,y+1]-phi[x,y-1])/2
	dxmz = (mz[x+1,y]-mz[x-1,y])/2
	dymz = (mz[x,y+1]-mz[x,y-1])/2

	l = ( - (mz[x,y] * (dxmz/mr) )*np.cos(phi[x,y])
		- mr*np.sin(phi[x,y])*dxphi
		- (mz[x,y] * (dymz/mr) )*np.sin(phi[x,y])
		+ mr*np.cos(phi[x,y])*dyphi )

	return l

def dphi_lagr(mz, phi, x, y):
	mr = np.sqrt(1 - (mz[x,y])**2)
	dxphi = (phi[x+1,y]-phi[x-1,y])/2
	dyphi = (phi[x,y+1]-phi[x,y-1])/2
	dxmz = (mz[x+1,y]-mz[x-1,y])/2
	dymz = (mz[x,y+1]-mz[x,y-1])/2

	dphi_l = ( (mz[x,y] * (dxmz/mr) )*np.sin(phi[x,y])
		- mr*np.cos(phi[x,y])*dxphi
		- (mz[x,y] * (dymz/mr) )*np.cos(phi[x,y])
		+ mr*np.sin(phi[x,y])*dyphi )

	return dphi_l

def dgradphi_lagr(mz, phi, x, y):
	mr = np.sqrt(1 - (mz[x,y])**2)

	dgradphi_l = np.array([-mr*np.sin(phi[x,y]), mr*np.cos(phi[x,y])])

	return dgradphi_l

def delta_phi_C(mz, phi, x, y):
	l = lagr(mz, phi, x, y)
	dphi_l = dphi_lagr(mz, phi, x, y)
	dgradphi_l = dgradphi_lagr(mz, phi, x, y)

	grad_ldgradphi_l =  ( ((lagr(mz, phi, x+1, y) * dgradphi_lagr(mz, phi, x+1, y)[0])
						-(lagr(mz, phi, x-1, y) * dgradphi_lagr(mz, phi, x-1, y)[0]))/2
						+((lagr(mz, phi, x, y+1) * dgradphi_lagr(mz, phi, x, y+1)[1])
						-(lagr(mz, phi, x, y-1) * dgradphi_lagr(mz, phi, x, y-1)[1]))/2 )

	return ( l*dphi_l - grad_ldgradphi_l )


def m_reconstruction_Dovzhenko(mz, num_steps, step_size):
	phi = np.zeros_like(mz)
	dlen = len(phi)
	tempphi = phi.copy()

	for i in range(0, num_steps):
		print(i)
		for x in range(2, dlen-2):
			for y in range(2, dlen-2):
				tempphi[x,y] = phi[x,y] - step_size*delta_phi_C(mz, phi, x, y)
		for x in range(2, dlen-2):
			for y in range(2, dlen-2):
				phi[x,y] =  tempphi[x,y]

	return phi
