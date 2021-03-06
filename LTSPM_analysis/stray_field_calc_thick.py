# @Author: Jenkins Alec <alec>
# @Date:   2017-01-18T09:54:01-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-09T15:04:19-07:00



import numpy as np
import fourier_image as fi

def stray_field_calc_thick(mx,my,mz,Ms,t,sim_size,z,windowData=False,windowPower=1/2):
    pi = np.pi
    mlen = len(mx)
    hlen = mlen/2
    res = sim_size/mlen
    fieldPrefactor = (4*pi*1.0e-7)*Ms

    # Mx = np.multiply(mx, 1/(res**2))
    if windowData:
        mx = fi.window_image(mx,windowPower)
        my = fi.window_image(my,windowPower)
        mz = fi.window_image(mz,windowPower)
    fmxs = np.fft.fft2(mx)
    fmys = np.fft.fft2(my)
    fmzs = np.fft.fft2(mz)

    fmx = np.fft.fftshift(fmxs)
    fmy = np.fft.fftshift(fmys)
    fmz = np.fft.fftshift(fmzs)

    kxv = 2*pi*np.linspace(-mlen/(2*sim_size), mlen/(2*sim_size), mlen, endpoint=False)
    kyv = 2*pi*np.linspace(-mlen/(2*sim_size), mlen/(2*sim_size), mlen, endpoint=False)
    kx, ky = np.meshgrid(kxv, kyv)

    k = np.sqrt(kx**2 + ky**2)

    alpha = np.divide(np.exp(-k*z),k+1e-10)*(1-np.exp(-k*t))/2

    Hzk = alpha*k*fmz - 1j*alpha*(kx*fmx + ky*fmy)

    surface_cdk = fmz
    volume_cdk = -1j*(kx*fmx + ky*fmy)

    Hxk = -1j * np.divide(kx, k+(1e-10)) * Hzk
    Hyk = -1j * np.divide(ky, k+(1e-10)) * Hzk

    meffk = surface_cdk+np.divide(volume_cdk, k+(1e-10))

    volume_cd = np.real(np.fft.ifft2(np.fft.ifftshift(volume_cdk)))
    surface_cd = np.real(np.fft.ifft2(np.fft.ifftshift(fmz)))
    meff = np.real(np.fft.ifft2(np.fft.ifftshift(meffk)))

    Hx = fieldPrefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hxk)))
    Hy = fieldPrefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hyk)))
    Hz = fieldPrefactor*np.real(np.fft.ifft2(np.fft.ifftshift(Hzk)))

    return [Hx, Hy, Hz], surface_cd, volume_cd, meff
