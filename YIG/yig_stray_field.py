# @Author: Jenkins Alec <alec>
# @Date:   2017-07-22T11:34:50-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-24T00:00:44-07:00



import numpy as np
import matplotlib.pyplot as plt
import json

import format_plots_tkagg as fp
import stray_field_calc_fast as sfcf
import stray_field_calc_thick as sfct

pi = np.pi

path = '/Users/alec/UCSB/LTSPM/YIG/'
material_params_path = path+'YIG_parameters.json'
with open(material_params_path, 'r') as fread:
    material_params = json.load(fread)

Ms = material_params['Ms']
t = material_params['t']

heights = np.array([5,10,20])*(1e-9)

r0s = np.array([5,10])*(1e-9)

dres = 0.5e-9
simSize = 200e-9
slen = int(np.round(simSize/dres))

x = np.linspace(-simSize, simSize, slen, endpoint=False)
y = np.linspace(-simSize, simSize, slen, endpoint=False)

xgrid, ygrid = np.meshgrid(x,y)

r = np.sqrt(xgrid**2 + ygrid**2)

bz_arrayB = np.zeros((2,3,slen,slen))
bz_arrayNL = np.zeros((2,3,slen,slen))
bz_arrayNR = np.zeros((2,3,slen,slen))

scd_arrayNR = np.zeros((2,3,slen,slen))
vcd_arrayNR = np.zeros((2,3,slen,slen))
scd_arrayB = np.zeros((2,3,slen,slen))
vcd_arrayB = np.zeros((2,3,slen,slen))
scd_arrayNL = np.zeros((2,3,slen,slen))
vcd_arrayNL = np.zeros((2,3,slen,slen))

for j in range(len(r0s)):

    # theta = 2*np.arctan(r0s[j]/(r + 1e-25))
    theta = pi * np.exp(-r/r0s[j]) / np.sqrt((2*r/r0s[j])+1)

    mz = np.cos(theta)
    mxB = np.sqrt(1-mz**2)*np.divide(-ygrid, np.sqrt(xgrid**2 + ygrid**2 + 1e-20))
    myB = np.sqrt(1-mz**2)*np.divide(xgrid, np.sqrt(xgrid**2 + ygrid**2 + 1e-20))
    mxNL = np.sqrt(1-mz**2)*np.divide(xgrid, np.sqrt(xgrid**2 + ygrid**2 + 1e-20))
    myNL = np.sqrt(1-mz**2)*np.divide(ygrid, np.sqrt(xgrid**2 + ygrid**2 + 1e-20))
    mxNR = np.sqrt(1-mz**2)*np.divide(-xgrid, np.sqrt(xgrid**2 + ygrid**2 + 1e-20))
    myNR = np.sqrt(1-mz**2)*np.divide(-ygrid, np.sqrt(xgrid**2 + ygrid**2 + 1e-20))

    for i in range(len(heights)):
        hNR, vcdNR, scdNR = sfct.stray_field_calc_thick(mxNR,myNR,mz,Ms,t,simSize,heights[i])
        hB, vcdB, scdB = sfct.stray_field_calc_thick(mxB,myB,mz,Ms,t,simSize,heights[i])
        hNL, vcdNL, scdNL = sfct.stray_field_calc_thick(mxNL,myNL,mz,Ms,t,simSize,heights[i])

        bz_arrayNR[j,i] = hNR[2]
        bz_arrayB[j,i] = hB[2]
        bz_arrayNL[j,i] = hNL[2]

        scd_arrayNR[j,i] = scdNR
        vcd_arrayNR[j,i] = vcdNR
        scd_arrayB[j,i] = scdB
        vcd_arrayB[j,i] = vcdB
        scd_arrayNL[j,i] = scdNL
        vcd_arrayNL[j,i] = vcdNL
# scd, vcd, meff, hk, h = sfcf.stray_field_calc_fast(mx,my,mz,Ms*t,simSize,heights[0])

cropSize = 100e-9
cropNum = int(slen*(cropSize/simSize))
cropRange = [int(slen/2-cropNum/2), int(slen/2+cropNum/2)]
cropExtent = (1e9)*cropSize/2
bzCropsB = bz_arrayB[:,:,cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]
bzCropsNL = bz_arrayNL[:,:,cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]
bzCropsNR = bz_arrayNR[:,:,cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]

scdCropsNR = scd_arrayNR[:,:,cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]
vcdCropsNR = vcd_arrayNR[:,:,cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]
scdCropsB = scd_arrayB[:,:,cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]
vcdCropsB = vcd_arrayB[:,:,cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]
scdCropsNL = scd_arrayNL[:,:,cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]
vcdCropsNL = vcd_arrayNL[:,:,cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]

mzCrop = mz[cropRange[0]:cropRange[1], cropRange[0]:cropRange[1]]

xCrop = np.linspace(-cropExtent, cropExtent, cropNum)

mphiB = np.arctan2(myB,mxB)

plt.close('all')

# fig1, ax1 = plt.subplots()
# plt.imshow(mxB, origin='lower', extent=[-cropSize/2, cropSize/2, -cropSize/2, cropSize/2])
# plt.colorbar()
#
# fig1, ax1 = plt.subplots()
# plt.imshow(myB, origin='lower', extent=[-cropSize/2, cropSize/2, -cropSize/2, cropSize/2])
# plt.colorbar()

# fig1, ax1 = plt.subplots()
# plt.imshow(mphiB, origin='lower', extent=[-cropSize/2, cropSize/2, -cropSize/2, cropSize/2])
# plt.colorbar()

fig1, ax1 = plt.subplots()
plt.imshow(bzCropsB[0,0], origin='lower', extent=[-cropExtent, cropExtent, -cropExtent, cropExtent])
plt.xlabel("x (nm)")
plt.ylabel("y (nm)")
plt.colorbar(fraction=0.046, pad=0.04, label=r"$B_z$ (T)")
plt.title(r"$B_z$, Bloch, $r_0$ = 5 nm, NV height = 5nm")

fig1, ax1 = plt.subplots()
plt.plot(xCrop, bzCropsB[0,0][:, int(cropNum/2)], label='NV height = 5nm')
plt.plot(xCrop, bzCropsB[0,1][:, int(cropNum/2)], label='NV height = 10nm')
plt.plot(xCrop, bzCropsB[0,2][:, int(cropNum/2)], label='NV height = 20nm')
plt.legend(loc=4,borderaxespad=1,prop={'size':10})
plt.xlabel("x (nm)")
plt.ylabel(r"$B_z$ (T)")
plt.title(r"$B_z$ linecuts, Bloch, $r_0$ = 5 nm")

fig1, ax1 = plt.subplots()
plt.plot(xCrop, bzCropsNR[0,1][:, int(cropNum/2)],  label=u'right-handed Néel')
plt.plot(xCrop, bzCropsB[0,1][:, int(cropNum/2)], label='Bloch')
plt.plot(xCrop, bzCropsNL[0,1][:, int(cropNum/2)], label=u'left-handed Néel')
plt.legend(loc=4,borderaxespad=1,prop={'size':9})
plt.xlabel("x (nm)")
plt.ylabel(r"$B_z$ (T)")
plt.title(r"$B_z$ linecuts, $r_0$ = 10 nm, NV height = 10 nm")

fig1, ax1 = plt.subplots()
plt.plot(xCrop, bzCropsB[0,1][:, int(cropNum/2)], label=r'$r_0$ = 5 nm')
plt.plot(xCrop, bzCropsB[1,1][:, int(cropNum/2)], label=r'$r_0$ = 10 nm')
plt.legend(loc=4,borderaxespad=1,prop={'size':10})
plt.xlabel("x (nm)")
plt.ylabel(r"$B_z$ (T)")
plt.title(r"$B_z$ linecuts, Bloch, NV height = 10nm")

fig1, ax1 = plt.subplots()
plt.plot(xCrop, scdCropsNR[1,1][:, int(cropNum/2)], label=r'$M_s \cdot \hat{z}$')
plt.plot(xCrop, (1e-8)*vcdCropsNR[1,1][:, int(cropNum/2)], label=r'($-\nabla \cdot M_s$) x thickness')
plt.legend(loc=1,borderaxespad=1,prop={'size':10})
plt.xlabel("x (nm)")
plt.ylabel("normalized effective surface charge density")
plt.title(r"$B_z$ linecuts, right-handed Néel, $r_0$ = 10 nm, NV height = 10nm")

fig1, ax1 = plt.subplots()
plt.plot(xCrop, scdCropsB[1,1][:, int(cropNum/2)], label=r'$M_s \cdot \hat{z}$')
plt.plot(xCrop, (1e-8)*vcdCropsB[1,1][:, int(cropNum/2)], label=r'($-\nabla \cdot M_s$) x thickness')
plt.legend(loc=4,borderaxespad=1,prop={'size':10})
plt.xlabel("x (nm)")
plt.ylabel("normalized effective surface charge density")
plt.title(r"$B_z$ linecuts, Bloch, $r_0$ = 10 nm, NV height = 10nm")

fig1, ax1 = plt.subplots()
plt.plot(xCrop, scdCropsNL[1,1][:, int(cropNum/2)], label=r'$M_s \cdot \hat{z}$')
plt.plot(xCrop, (1e-8)*vcdCropsNL[1,1][:, int(cropNum/2)], label=r'($-\nabla \cdot M_s$) x thickness')
plt.legend(loc=4,borderaxespad=1,prop={'size':10})
plt.xlabel("x (nm)")
plt.ylabel("normalized effective surface charge density")
plt.title(r"$B_z$ linecuts, left-handed Néel, $r_0$ = 10 nm, NV height = 10nm")

# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, bzCropsNL[0,0][:, int(cropNum/2)], label='5 nm NV height')
# plt.plot(xCrop, bzCropsNL[0,1][:, int(cropNum/2)], label='10 nm NV height')
# plt.plot(xCrop, bzCropsNL[0,2][:, int(cropNum/2)], label='20 nm NV height')
# plt.legend(loc=3,borderaxespad=1,prop={'size':10})
# plt.title(u"Bz linecuts, left-handed Néel, 5 nm r0")

# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, bzCropsNR[0,0][:, int(cropNum/2)], label='5 nm NV height')
# plt.plot(xCrop, bzCropsNR[0,1][:, int(cropNum/2)], label='10 nm NV height')
# plt.plot(xCrop, bzCropsNR[0,2][:, int(cropNum/2)], label='20 nm NV height')
# plt.legend(loc=3,borderaxespad=1,prop={'size':10})
# plt.title(u"Bz linecuts, right-handed Néel, 5 nm r0")
#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, bzCropsNR[0,0][:, int(cropNum/2)],  label=u'right-handed Néel')
# plt.plot(xCrop, bzCropsB[0,1][:, int(cropNum/2)], label='Bloch')
# plt.plot(xCrop, bzCropsNL[0,2][:, int(cropNum/2)], label=u'left-handed Néel')
# plt.legend(loc=3,borderaxespad=1,prop={'size':10})
# plt.title("Bz linecuts, 5 nm r0, 10 nm NV height")

fp.format_plots(plt, small=True)

plt.show()
