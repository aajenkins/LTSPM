# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T09:54:01-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-10T13:53:32-07:00



import numpy as np
import fourier_image as fi

def stray_field_calc_fast(mx,my,mz,Mst,sim_size,z,windowData=False):
    pi = np.pi
    mlen = len(mx)
    hlen = mlen/2
    res = sim_size/mlen
    fieldPrefactor = (4*pi*1.0e-7)*Mst/2

    # Mx = np.multiply(mx, 1/(res**2))
    if windowData:
        mx = fi.window_image(mx,windowPower)
        my = fi.window_image(my,windowPower)
        mz = fi.window_image(mz,windowPower)
    fmxs = np.fft.fft2(mx, norm="ortho")
    fmys = np.fft.fft2(my, norm="ortho")
    fmzs = np.fft.fft2(mz, norm="ortho")

    fmx = np.fft.fftshift(fmxs)
    fmy = np.fft.fftshift(fmys)
    fmz = np.fft.fftshift(fmzs)

    surface_cdk = fmz

    volume_cdk = np.zeros_like(fmx)
    meffk = np.zeros_like(fmx)

    kxv = np.arange(-pi*mlen/sim_size, pi*mlen/sim_size, 2*pi/sim_size)
    kyv = np.arange(-pi*mlen/sim_size, pi*mlen/sim_size, 2*pi/sim_size)
    kx, ky = np.meshgrid(kxv, kyv)

    k = np.sqrt(kx**2 + ky**2)

    volume_cdk = -1j*(kx*fmx + ky*fmy)

    Hxk = np.zeros_like(fmx)
    Hyk = np.zeros_like(fmx)
    Hzk = np.zeros_like(fmx)

    Hzk = np.exp(-k*z) * k * (surface_cdk + np.divide(volume_cdk, k))
    Hxk = -1j * np.divide(kx, k) * Hzk
    Hyk = -1j * np.divide(ky, k) * Hzk
    meffk = surface_cdk+np.divide(volume_cdk, k)

    surface_cd = np.real(np.fft.ifft2(np.fft.ifftshift(surface_cdk), norm="ortho"))
    volume_cd = np.real(np.fft.ifft2(np.fft.ifftshift(volume_cdk), norm="ortho"))
    meff = np.real(np.fft.ifft2(np.fft.ifftshift(meffk), norm="ortho"))

    Hx = fieldPrefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hxk), norm="ortho"))
    Hy = fieldPrefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hyk), norm="ortho"))
    Hz = fieldPrefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hzk), norm="ortho"))

    return surface_cd, volume_cd, meff, [Hxk, Hyk, Hzk], [Hx, Hy, Hz]
