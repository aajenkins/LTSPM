# @Author: Jenkins Alec <alec>
# @Date:   2017-02-03T18:47:38-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-04T19:38:58-08:00



import numpy as np
from scipy import signal

pi = np.pi

def f(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    fp1, fp2, fp3, fp4 = 0., 0., 0., 0.
    fp1 = (r/6.)*(3.*(z**2) - r**2)
    if (y != 0 and z != 0):
        fp2 = 0.5 * ( x * (y**2 - z**2)*np.arcsinh(x/np.sqrt(y**2 + z**2)) )
    if (x != 0 and z != 0):
        fp3 = 0.5 * ( y * (x**2 - z**2)*np.arcsinh(y/np.sqrt(x**2 + z**2)) )
    if (z>0):
        fp4 = -x*y*z*np.arctan(x*y/(z*r))
    f = fp1 + fp2 + fp3 + fp4
    return f

def F2(x,y,z):
    return ( f(x,y,z) - f(0,y,z) - f(x,0,z) + f(0,0,z) )

def F1(x,y,z,a):
    return ( F2(x,y,z) - F2(x,y-a,z)
            - F2(x-a,y,z) + F2(x-a,y-a,z) )

def F(x,y,z,a):
    return ( F1(x+a,y+a,z,a) - F1(x+a,y,z,a)
            - F1(x,y+a,z,a) + F1(x,y,z,a) )

def Nzz(m,n,t,a):
    if (m == 0 and n == 0):
        Nzzmn = 0
    else:
        Nzzmn = (1/(4*pi*a*a*t)) * ( 2*F(m*a,n*a,0,a)
                                    - F(m*a,n*a,t,a) - F(m*a,n*a,-t,a) )
    return Nzzmn

def calc_demag_tensor_newell(slen, a, t, path='demag_tensor.txt'):
    d11 = np.zeros((slen,slen))
    for j in range(0, slen):
        for i in range(0, slen):
            d11[j,i] = Nzz(i,j,t,a)

    d10 = np.flipud(d11[1:,:])
    d00 = np.fliplr(np.flipud(d11[1:,1:]))
    d01 = np.fliplr(d11[:,1:])

    d = np.concatenate((np.concatenate((d00,d10), axis=1),
                       np.concatenate((d01,d11), axis=1)), axis=0)

    np.savetxt(path, d, delimiter=',')
