# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T09:54:01-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-23T19:25:17-07:00



import numpy as np
import fourier_image as fi

def stray_field_calc_thick(mx,my,mz,Ms,t,sim_size,z,windowData=False, windowPower=1/2):
    pi = np.pi
    mlen = len(mx)
    hlen = mlen/2
    res = sim_size/mlen
    fieldPrefactor = (4*pi*1.0e-7)*Ms/2

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

    kxv = np.arange(-pi*mlen/sim_size, pi*mlen/sim_size, 2*pi/sim_size)
    kyv = np.arange(-pi*mlen/sim_size, pi*mlen/sim_size, 2*pi/sim_size)
    kx, ky = np.meshgrid(kxv, kyv)

    k = np.sqrt(kx**2 + ky**2)

    # alphazk = np.divide(np.exp(-(z+t)*k) * (np.exp(t*k)-1), k**2)
    # alphaxyk = np.divide(np.exp(-(z+t)*k) * (np.exp(t*k)-1), k)
    alphazk = np.divide(t*np.exp(-z*k), k)
    alphaxyk = t*np.exp(-z*k)

    Hzk = alphazk*(k**2)*fmz - 1j*alphaxyk*(kx*fmx + ky*fmy)

    surface_cdk = fmz
    volume_cdk = -1j*(kx*fmx + ky*fmy)

    Hxk = -1j * np.divide(kx, k) * Hzk
    Hyk = -1j * np.divide(ky, k) * Hzk

    volume_cd = np.real(np.fft.ifft2(np.fft.ifftshift(volume_cdk), norm="ortho"))
    surface_cd = np.real(np.fft.ifft2(np.fft.ifftshift(surface_cdk), norm="ortho"))

    Hx = fieldPrefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hxk), norm="ortho"))
    Hy = fieldPrefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hyk), norm="ortho"))
    Hz = fieldPrefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hzk), norm="ortho"))

    return [Hx, Hy, Hz], volume_cd, surface_cd
