# @Author: Jenkins Alec <alec>
# @Date:   2017-01-20T18:13:22-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-09T14:09:36-07:00



import numpy as np
import matplotlib.pyplot as plt
import stray_field_calc_thick as sfct
import plotting.format_plots_tkagg as fp
import json

import fourier_image as fi

pi = np.pi
scannum = 1760

path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
material_params_path = path+'material_parameters.json'
with open(material_params_path, 'r') as fread:
    material_params = json.load(fread)

Ms = material_params['Ms']
MstError = material_params['MstError']
t = material_params['t']

mu0 = 4*pi*(1e-7)

B0 = (1e4)*mu0*Ms*t/(2*pi)
height = 1e-7
size = 2.5e-4

def edge_field(x, *args):
    h = args[0]
    x0 = 0
    bx = B0*(h/(h**2+(x-x0)**2))
    bz = -B0*((x-x0)/(h**2+(x-x0)**2))

    return bx, bz

baspath = "/Users/alec/UCSB/cofeb_analysis_data/stray_field_test/"
filenames = ["medgex_test","medgey_test","medgez_test"]
numfiles = len(filenames)
m = []
for i in range(0,numfiles):
    m.append(np.transpose(np.loadtxt(baspath+filenames[i]+".dat")))

mwx = fi.window_image(m[1])
mwy = fi.window_image(m[0])
mwz = fi.window_image(m[2])
h,_,_,_ = sfct.stray_field_calc_thick(mwx,mwy,mwz,Ms,t,size,height)

h[0]=(1e4)*h[0]
h[1]=(1e4)*h[1]
h[2]=(1e4)*h[2]

flen = len(h[0])

fx = np.linspace(-size/2, size/2, flen, endpoint=False)


bx, bz = edge_field(fx, height)

plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(fx, h[1][:,int(flen/2)])
plt.plot(fx, bx)

fig1, ax1 = plt.subplots()
plt.plot(fx, h[2][:,int(flen/2)])
plt.plot(fx, bz)

fp.format_plots(plt)

plt.show()
