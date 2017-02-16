# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T12:05:35-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-12T15:12:34-08:00

import numpy as np

pi = np.pi

def d_phi(phi1, phi0):
	if ((phi1-phi0)%(2*pi) > pi):
	    dphi = (phi1-phi0)%(2*pi)-2*pi
	else:
	    dphi = (phi1-phi0)%(2*pi)
	return dphi

def lagr(mz, phi, x, y):
	mr = np.sqrt(1 - (mz[x,y])**2)
	dxphi = d_phi(phi[x+1,y],phi[x-1,y])/2
	dyphi = d_phi(phi[x,y+1],phi[x,y-1])/2
	dxmz = (mz[x+1,y]-mz[x-1,y])/2
	dymz = (mz[x,y+1]-mz[x,y-1])/2

	l = ( -(mz[x,y]/mr)*((dxmz)*np.cos(phi[x,y]) + (dymz)*np.sin(phi[x,y]))
		+ mr*(-dxphi*np.sin(phi[x,y]) + dyphi*np.cos(phi[x,y])) )

	return l

def dphi_lagr(mz, phi, x, y):
	mr = np.sqrt(1 - (mz[x,y])**2)
	dxphi = d_phi(phi[x+1,y],phi[x-1,y])/2
	dyphi = d_phi(phi[x,y+1],phi[x,y-1])/2
	dxmz = (mz[x+1,y]-mz[x-1,y])/2
	dymz = (mz[x,y+1]-mz[x,y-1])/2

	dphi_l = ( -(mz[x,y]/mr)*(-(dxmz)*np.sin(phi[x,y]) + (dymz)*np.cos(phi[x,y]))
		+ mr*(-dxphi*np.cos(phi[x,y]) - dyphi*np.sin(phi[x,y])) )

	return dphi_l

def dgradphi_lagr(mz, phi, x, y):
	mr = np.sqrt(1 - (mz[x,y])**2)

	dgradphi_l = np.array([-mr*np.sin(phi[x,y]), mr*np.cos(phi[x,y])])

	return dgradphi_l

def delta_phi_C(mz, phi, x, y):
	l = lagr(mz, phi, x, y)
	dphi_l = dphi_lagr(mz, phi, x, y)
	grad_lx = (lagr(mz, phi, x+1, y)-lagr(mz, phi, x-1, y))/2
	grad_ly = (lagr(mz, phi, x, y+1)-lagr(mz, phi, x, y-1))/2
	grad_l_dot_dgradphi_lagr = ( grad_lx * dgradphi_lagr(mz, phi, x, y)[0]
								+ grad_ly * dgradphi_lagr(mz, phi, x, y)[1] )
	# grad_ldgradphi_l =  ( ((lagr(mz, phi, x+1, y) * dgradphi_lagr(mz, phi, x+1, y)[0])
	# 					-(lagr(mz, phi, x-1, y) * dgradphi_lagr(mz, phi, x-1, y)[0]))/2
	# 					+((lagr(mz, phi, x, y+1) * dgradphi_lagr(mz, phi, x, y+1)[1])
	# 					-(lagr(mz, phi, x, y-1) * dgradphi_lagr(mz, phi, x, y-1)[1]))/2 )
	return ( l*dphi_l - grad_l_dot_dgradphi_lagr )


def m_reconstruction_Dovzhenko(mz, phi_seed, max_num_steps, cutoff_residual, max_angle_step):
	mzlen = len(mz)
	phi = phi_seed.copy()
	# phileft = np.full((mzlen, int(mzlen/2)), np.pi)
	# phi[:, 0:int(mzlen/2)] = phileft
	dlen = len(phi)
	phi_gradient = np.zeros_like(phi)
	step_size = max_angle_step
	i = 0
	total_residual = cutoff_residual+1
	max_abs_gradient = cutoff_residual+1

	while (max_abs_gradient > cutoff_residual and i < max_num_steps):
		for x in range(2, dlen-2):
			for y in range(2, dlen-2):
				phi_gradient[x,y] = delta_phi_C(mz, phi, x, y)
		max_gradient = np.max(phi_gradient)
		min_gradient = np.min(phi_gradient)
		max_abs_gradient = np.max([max_gradient, np.abs(min_gradient)])
		step_size = max_angle_step/max_abs_gradient
		total_residual = np.sum(phi_gradient**2)

		for x in range(2, dlen-2):
			for y in range(2, dlen-2):
				phi[x,y] = phi[x,y] - step_size*phi_gradient[x,y]

		if (i%5 == 0):
			print(str(i)+', residual = {:04f}, max/min = {:04f} / {:04f}'.format(total_residual, max_gradient, min_gradient))

		i = i+1

	return phi
