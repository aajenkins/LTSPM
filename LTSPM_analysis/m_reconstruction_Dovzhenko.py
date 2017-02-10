# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T12:05:35-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-10T14:44:24-08:00



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


def m_reconstruction_Dovzhenko(mz, max_num_steps, cutoff_residual, max_angle_step):
	phi = np.zeros_like(mz)
	dlen = len(phi)
	tempphi = phi.copy()
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

		print(str(i)+', residual = '+
		str(total_residual)+', max/min = '+str(max_gradient)+' / '
		+str(min_gradient))

		i = i+1

	return phi
