# @Author: Jenkins Alec <alec>
# @Date:   2017-02-03T18:47:38-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-04T15:45:43-08:00



import numpy as np
from scipy import signal

pi = np.pi

def f220(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    fp1 = 0
    fp2 = 0
    fp3 = 0
    fp1 = (r/6)*(3*(z**2) - r**2)
    if (r != 0.0 and np.abs(x) != np.abs(r) and np.abs(y) != np.abs(r)):
        fp2 = 0.5 * ( x * (y**2 - z**2)*np.arctanh(x/r)
                    + y * (x**2 - z**2)*np.arctanh(y/r) )
    if (z>0):
        fp3 = -x*y*z*np.arctan(x*y/(z*r))
    f = fp1 + fp2 + fp3
    return f

def frr(m,n,t,a):
    return ( 2 * (f220(m*a,n*a,t) - f220(m*a,n*a,0)) )

def frrp(m,n,t,a):
    return ( frr(m,n+1,t,a) + frr(m,n-1,t,a) - 2*frr(m,n,t,a) )

def d_interaction(m,n,a,t):
    if (m == 0 and n == 0):
        dmn = 0
    else:
        dmn = (1/(4*pi*a*a*t))*( frrp(m+1,n,t,a) + frrp(m-1,n,t,a)
                       - 2*frrp(m,n,t,a) )
    return dmn

def calc_demag_tensor(slen, a, t, path='demag_tensor.txt'):
    d11 = np.zeros((slen,slen))
    for j in range(0, slen):
        for i in range(0, slen):
            d11[j,i] = d_interaction(i,j,a,t)

    d10 = np.flipud(d11[1:,:])
    d00 = np.fliplr(np.flipud(d11[1:,1:]))
    d01 = np.fliplr(d11[:,1:])

    d = np.concatenate((np.concatenate((d00,d10), axis=1),
                       np.concatenate((d01,d11), axis=1)), axis=0)

    np.savetxt(path, d, delimiter=',')
