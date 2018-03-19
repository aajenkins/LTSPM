# @Author: Jenkins Alec <alec>
# @Date:   2018-03-13T22:17:55-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2018-03-13T23:14:15-07:00



from scipy.optimize import curve_fit
from scipy import signal
import numpy as np

def j_gurzhi(x, i0, w, l):
	j = (i0/(w-2*l*np.tanh(w/(2*l)))) * ( 1- (np.cosh(x/l)/np.cosh(w/(2*l))) )
	return j

def current_strip_field_bx(x, x0, z, i0, w, l):
    bx = (2e-7) * j_gurzhi(x, i0, w, l) * (z / (z**2 + (x0-x)**2) )
    return bx

def current_strip_field_bz(x, x0, z, i0, w, l):
    bz = (2e-7) * j_gurzhi(x, i0, w, l) * ((x0-x) / (z**2 + (x0-x)**2) )
    return bz
