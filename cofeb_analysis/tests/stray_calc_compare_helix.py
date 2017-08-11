# @Author: Jenkins Alec <alec>
# @Date:   2017-01-20T18:13:22-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-09T14:57:20-07:00



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
t = 1e-8

mu0 = 4*pi*(1e-7)

B0 = mu0*Ms/2
height = 5e-8
l=5e-8
res = 200
size = l*10

def helix_field(x, *args):
    h = args[0]
    bx = B0*(np.exp(-2*pi*h/l)-np.exp(-2*pi*(h+t)/l))*np.sin(2*pi*x/l)
    bz = B0*(np.exp(-2*pi*h/l)-np.exp(-2*pi*(h+t)/l))*np.cos(2*pi*x/l)

    return bx, bz


x = np.linspace(-size/2, size/2, res, endpoint=False)
y = x
xx,yy = np.meshgrid(x,y)

my = np.sin(2*pi*xx/l)
mz = np.cos(2*pi*xx/l)
mx = np.zeros_like(my)

mwx = fi.window_image(mx, power=1/20)
mwy = fi.window_image(my, power=1/20)
mwz = fi.window_image(mz, power=1/20)
h,_,_,_ = sfct.stray_field_calc_thick(mx,my,mz,Ms,t,size,height)

flen = len(h[0])

fx = np.linspace(-size/2, size/2, flen, endpoint=False)

bx, bz = helix_field(fx, height)

plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(fx, h[0][int(flen/2),:])
plt.plot(fx, bx)

fig1, ax1 = plt.subplots()
plt.plot(fx, h[2][int(flen/2),:])
plt.plot(fx, bz)

fp.format_plots(plt)

plt.show()
