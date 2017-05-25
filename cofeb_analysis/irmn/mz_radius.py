# @Author: Jenkins Alec <alec>
# @Date:   2017-03-26T15:54:40-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-03-26T16:04:40-07:00



import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import misc
from scipy import signal
from scipy.optimize import curve_fit
import json
import matplotlib.pylab as pylab
import vector_reconstruction as vr
import load_scan as lscan
import fourier_image as fi
import format_plot as fp

pi=np.pi
path = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
mzdata = np.loadtxt(path+'mz.dat')
scanL = 2.5e-6
xres = 50

def fit_tanh(x, *params):
    y = np.zeros_like(x)
    c = params[0]
    a = params[1]
    x0 = params[2]
    wid = params[3]

    y = c+(a/2)*np.tanh((x-x0)/wid)
    return y

x0, y0 = 22, 28
phinum = 16
lcnum = 15
lclen = 15
mphi = np.zeros((phinum,lcnum))
for i in range(0,phinum):
    phi = i*2*pi/phinum
    x1, y1 = x0-lclen*np.cos(phi), y0-lclen*np.sin(phi)
    x, y = np.linspace(x0, x1, lcnum), np.linspace(y0, y1, lcnum)
    mphi[i] = ndimage.map_coordinates(np.transpose(mzdata), np.vstack((x,y)), order=1)

xf = np.arange(0,lcnum)
fits = np.zeros((phinum,lcnum))
guesses = np.zeros((phinum, 4))
widths = np.zeros(phinum)
r0s = np.zeros(phinum)
angles = np.linspace(0,2*pi,phinum)

for i in range (0,phinum):
	y = mphi[i]
	guesses[i] = [(y[-1]+y[0])/2,y[-1]-y[0],6,1]
	popt, pcov = curve_fit(fit_tanh, xf, mphi[i], p0=guesses[i])
	fits[i] = fit_tanh(xf, *popt)
	widths[i] = np.abs(popt[3])
	r0s[i] = popt[2]*scanL/xres

np.savetxt('/Users/alec/UCSB/cofeb_analysis_data/irmn/radiusphi.txt',(angles,r0s), delimiter=',')

plt.close('all')

fig, ax = plt.subplots()
im = plt.imshow(mzdata)
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

fig, axes = plt.subplots(nrows=phinum, sharex=True, sharey=True)
fig.set_size_inches(3, 5)

for i in range(0,phinum):
	axes[i].plot(mphi[i],'b.')
	axes[i].plot(xf, fit_tanh(xf, *guesses[i]), 'g')
	axes[i].plot(xf, fits[i], 'r')
	axes[i].get_yaxis().set_visible(False)

fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
#fig.tight_layout()
# ax1.xaxis.set_ticklabels([])
# ax1.yaxis.set_ticklabels([])
fp.format_plot(plt, 400, 900, 900, 50, tight=False)

plt.show()
