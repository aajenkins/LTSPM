# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T09:54:01-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-16T13:11:07-08:00



import numpy as np
import fourier_image as fi
# from scipy import fftpack

def stray_field_calc_discrete(mx,my,mz,Mst,sim_size,z):
    pi = np.pi
    mlen = len(mx)
    hlen = mlen/2
    res = sim_size/mlen
    sprefactor = Mst

    # Mx = np.multiply(mx, 1/(res**2))
    rpotential = np.zeros_like(mx)
    rpotential_k = np.zeros_like(mx)

    fmxs = np.fft.fft2(mx, norm="ortho")
    fmys = np.fft.fft2(my, norm="ortho")
    fmzs = np.fft.fft2(mz, norm="ortho")

    fmx = np.fft.fftshift(fmxs)
    fmy = np.fft.fftshift(fmys)
    fmz = np.fft.fftshift(fmzs)

    surface_cdk = fmz

    volume_cdk = np.zeros_like(fmx)
    meffk = np.zeros_like(fmx)

    for j in range(0,mlen):
        ky = 2*pi*(j-hlen)/sim_size
        y = res*(j-hlen)
        for i in range(0, mlen):
            kx = 2*pi*(i-hlen)/sim_size
            x = res*(i-hlen)
            rpotential[j][i] = 1/np.sqrt(x**2 + y**2 + z**2)
            volume_cdk[j][i] = -1j*(kx*fmx[j,i]+ky*fmy[j,i])


    Hxk = np.zeros_like(fmx)
    Hyk = np.zeros_like(fmx)
    Hzk = np.zeros_like(fmx)

    rpotential_k = np.fft.fftshift(np.fft.fft2(rpotential, norm="ortho"))

    for j in range(0, mlen):
        ky = 2 * pi * (j-hlen) / sim_size
        for i in range(0, mlen):
            kx = 2 * pi * (i-hlen) / sim_size
            k = np.sqrt(kx**2 + ky**2)
            if (k==0):
                Hzk[j,i] = 0
                Hxk[j][i] = 0
                Hyk[j][i] = 0
            else:
                Hzk[j][i] = (1/(4*pi))*rpotential[j][i]*(surface_cdk[j,i]+(volume_cdk[j][i])/k)
                Hxk[j][i] = -1j * (kx / k) * Hzk[j, i]
                Hyk[j][i] = -1j * (ky / k) * Hzk[j, i]
                meffk[j][i] = surface_cdk[j,i]+(volume_cdk[j][i])/k

    surface_cd = np.real(np.fft.ifft2(np.fft.ifftshift(surface_cdk), norm="ortho"))
    volume_cd = np.real(np.fft.ifft2(np.fft.ifftshift(volume_cdk), norm="ortho"))
    meff = np.real(np.fft.ifft2(np.fft.ifftshift(meffk), norm="ortho"))

    Hx = sprefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hxk), norm="ortho"))
    Hy = sprefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hyk), norm="ortho"))
    Hz = sprefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hzk), norm="ortho"))

    return surface_cd, volume_cd, meff, [Hxk, Hyk, Hzk], [Hx, Hy, Hz]
