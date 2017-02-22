# @Author: Jenkins Alec <alec>
# @Date:   2017-02-16T12:00:08-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-16T12:40:59-08:00


import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
res = 3
z = 5
slen = 100
hlen = int(slen/2)

f = np.zeros((slen, slen))
fkcont = np.zeros((slen, slen))

for j in range(slen):
    y = res * (j-hlen)
    kj = 2*pi*(j-hlen)/(res*slen)
    for i in range(slen):
        x = res * (i-hlen)
        ki = 2*pi*(i-hlen)/(res*slen)
        k = np.sqrt(ki**2 + kj**2)
        f[j][i] = (1/np.sqrt(x**2 + y**2 + z**2))
        if (k != 0):
            fkcont[j][i] = (1/k)*np.exp(-k*z)

fk = np.fft.fftshift(np.fft.fft2(f))
frecon = np.real(np.fft.ifft2(fk, norm='ortho'))
freconcont = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.fftshift(fkcont), norm='ortho')))

plt.close('all')

fig, ax = plt.subplots()
im = plt.imshow(np.abs(fk))
plt.colorbar()

plt.show()
