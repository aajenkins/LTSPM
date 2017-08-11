# @Author: Jenkins Alec <alec>
# @Date:   2017-07-22T11:34:50-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-03T12:44:09-07:00



import numpy as np
import matplotlib.pyplot as plt
import json

import plotting.format_plots_tkagg as fp
import stray_field_calc_fast as sfcf
import stray_field_calc_thick as sfct

pi = np.pi

path = '/Users/alec/UCSB/LTSPM/YIG/'
material_params_path = path+'YIG_parameters.json'
with open(material_params_path, 'r') as fread:
    material_params = json.load(fread)

Ms = material_params['Ms']
t = material_params['t']

heights = np.array([10])*(1e-9)
# heights = np.array([5,10,20,30,40,50])*(1e-9)
heightlen = len(heights)

FWHM = np.array([5,10,20])*(1e-9)
r0s = 2*FWHM

r0slen = len(r0s)

dres = 1.0e-9
simSize = 200e-9
slen = int(np.round(2*simSize/dres))

x = np.linspace(-simSize, simSize, slen, endpoint=False)
y = np.linspace(-simSize, simSize, slen, endpoint=False)

xgrid, ygrid = np.meshgrid(x,y)

r = np.sqrt(xgrid**2 + ygrid**2)

bz_arrayB = np.zeros((r0slen,heightlen,slen,slen))
bz_arrayNL = np.zeros((r0slen,heightlen,slen,slen))
bz_arrayNR = np.zeros((r0slen,heightlen,slen,slen))

scd_arrayNR = np.zeros((r0slen,heightlen,slen,slen))
vcd_arrayNR = np.zeros((r0slen,heightlen,slen,slen))
scd_arrayB = np.zeros((r0slen,heightlen,slen,slen))
vcd_arrayB = np.zeros((r0slen,heightlen,slen,slen))
scd_arrayNL = np.zeros((r0slen,heightlen,slen,slen))
vcd_arrayNL = np.zeros((r0slen,heightlen,slen,slen))

fbz_arrayB = np.zeros((r0slen,heightlen,slen,slen), dtype=np.complex)

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
        hNR, scdNR, vcdNR, _ = sfct.stray_field_calc_thick(mxNR,myNR,mz,Ms,t,simSize,heights[i])
        hB, scdB, vcdB, fhB = sfct.stray_field_calc_thick(mxB,myB,mz,Ms,t,simSize,heights[i])
        hNL, scdNL, vcdNL, _ = sfct.stray_field_calc_thick(mxNL,myNL,mz,Ms,t,simSize,heights[i])

        bz_arrayNR[j,i] = hNR[2]
        bz_arrayB[j,i] = hB[2]
        bz_arrayNL[j,i] = hNL[2]

        scd_arrayNR[j,i] = scdNR
        vcd_arrayNR[j,i] = vcdNR
        scd_arrayB[j,i] = scdB
        vcd_arrayB[j,i] = vcdB
        scd_arrayNL[j,i] = scdNL
        vcd_arrayNL[j,i] = vcdNL

        fbz_arrayB[j,i] = fhB
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

savepath = "/Users/alec/UCSB/updates/YIG_7_23_17/"

plt.close('all')

# fig1, ax1 = plt.subplots()
# plt.imshow(bzCropsB[0,0], origin='lower', extent=[-cropExtent, cropExtent, -cropExtent, cropExtent])
# plt.xlabel("x (nm)")
# plt.ylabel("y (nm)")
# plt.colorbar(fraction=0.046, pad=0.04, label=r"$B_z$ (T)")
# plt.title(r"$B_z$, Bloch, $r_0$ = 5 nm, NV height = 5nm")
# plt.savefig(savepath+'bloch_image_5nm_r05nm.pdf',  bbox_inches='tight')
#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, bzCropsB[0,0][:, int(cropNum/2)], label='NV height = 5nm')
# plt.plot(xCrop, bzCropsB[0,1][:, int(cropNum/2)], label='NV height = 10nm')
# plt.plot(xCrop, bzCropsB[0,2][:, int(cropNum/2)], label='NV height = 20nm')
# plt.legend(loc=4,borderaxespad=1,prop={'size':10})
# plt.xlabel("x (nm)")
# plt.ylabel(r"$B_z$ (T)")
# plt.title(r"$B_z$ linecuts, Bloch, $r_0$ = 5 nm")
# plt.savefig(savepath+'bloch_cut_heights_r05nm.pdf',  bbox_inches='tight')
#
#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, bzCropsNR[1,1][:, int(cropNum/2)],  label=u'right-handed Néel')
# plt.plot(xCrop, bzCropsB[1,1][:, int(cropNum/2)], label='Bloch')
# plt.plot(xCrop, bzCropsNL[1,1][:, int(cropNum/2)], label=u'left-handed Néel')
# plt.legend(loc=4,borderaxespad=1,prop={'size':9})
# plt.xlabel("x (nm)")
# plt.ylabel(r"$B_z$ (T)")
# plt.title(r"$B_z$ linecuts, $r_0$ = 10 nm, NV height = 10 nm")
# plt.savefig(savepath+'dw_comp_10nm_r010nm.pdf',  bbox_inches='tight')
#
fig1, ax1 = plt.subplots()
plt.plot(xCrop, bzCropsNR[0,0][:, int(cropNum/2)], label=r'FWHM = 5 nm')
plt.plot(xCrop, bzCropsNR[1,0][:, int(cropNum/2)], label=r'FWHM = 10 nm')
plt.plot(xCrop, bzCropsNR[2,0][:, int(cropNum/2)], label=r'FWHM = 20 nm')
plt.legend(loc=1,borderaxespad=1,prop={'size':12})
plt.xlabel("x (nm)")
plt.ylabel(r"$B_z$ (T)")
plt.title(r"$B_z$ linecuts, right-handed Néel, 10nm NV height")
plt.savefig(savepath+'rn_fwhm_10nm.pdf',  bbox_inches='tight')

fig1, ax1 = plt.subplots()
plt.plot(xCrop, bzCropsNL[0,0][:, int(cropNum/2)], label=r'FWHM = 5 nm')
plt.plot(xCrop, bzCropsNL[1,0][:, int(cropNum/2)], label=r'FWHM = 10 nm')
plt.plot(xCrop, bzCropsNL[2,0][:, int(cropNum/2)], label=r'FWHM = 20 nm')
plt.legend(loc=4,borderaxespad=1,prop={'size':12})
plt.xlabel("x (nm)")
plt.ylabel(r"$B_z$ (T)")
plt.title(r"$B_z$ linecuts, left-handed Néel, 10nm NV height")
plt.savefig(savepath+'ln_fwhm_10nm.pdf',  bbox_inches='tight')
#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, scdCropsNR[1,1][:, int(cropNum/2)], label=r'$M_s \cdot \hat{z}$')
# plt.plot(xCrop, (1e-8)*vcdCropsNR[1,1][:, int(cropNum/2)], label=r'($-\nabla \cdot M_s$) x thickness')
# plt.legend(loc=1,borderaxespad=1,prop={'size':10})
# plt.xlabel("x (nm)")
# plt.ylabel("normalized effective surface charge density")
# plt.title(r"$B_z$ linecuts, right-handed Néel, $r_0$ = 10 nm, NV height = 10nm")
# plt.savefig(savepath+'charge_density_NR.pdf',  bbox_inches='tight')
#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, scdCropsB[1,1][:, int(cropNum/2)], label=r'$M_s \cdot \hat{z}$')
# plt.plot(xCrop, (1e-8)*vcdCropsB[1,1][:, int(cropNum/2)], label=r'($-\nabla \cdot M_s$) x thickness')
# plt.legend(loc=4,borderaxespad=1,prop={'size':10})
# plt.xlabel("x (nm)")
# plt.ylabel("normalized effective surface charge density")
# plt.title(r"$B_z$ linecuts, Bloch, FWHM = 10 nm, NV height = 10nm")
# plt.savefig(savepath+'charge_density_B.pdf',  bbox_inches='tight')
#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, scdCropsNL[1,1][:, int(cropNum/2)], label=r'$M_s \cdot \hat{z}$')
# plt.plot(xCrop, (1e-8)*vcdCropsNL[1,1][:, int(cropNum/2)], label=r'($-\nabla \cdot M_s$) x thickness')
# plt.legend(loc=4,borderaxespad=1,prop={'size':10})
# plt.xlabel("x (nm)")
# plt.ylabel("normalized effective surface charge density")
# plt.title(r"$B_z$ linecuts, left-handed Néel, $r_0$ = 10 nm, NV height = 10nm")
# plt.savefig(savepath+'charge_density_NL.pdf',  bbox_inches='tight')
#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, bzCropsNL[0,0][:, int(cropNum/2)], label='5 nm NV height')
# plt.plot(xCrop, bzCropsNL[0,1][:, int(cropNum/2)], label='10 nm NV height')
# plt.plot(xCrop, bzCropsNL[0,2][:, int(cropNum/2)], label='20 nm NV height')
# plt.legend(loc=3,borderaxespad=1,prop={'size':10})
# plt.title(u"Bz linecuts, left-handed Néel, 5 nm r0")

# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, bzCropsNR[1,0][:, int(cropNum/2)], label='5 nm NV height')
# plt.plot(xCrop, bzCropsNR[1,1][:, int(cropNum/2)], label='10 nm NV height')
# plt.plot(xCrop, bzCropsNR[1,2][:, int(cropNum/2)], label='20 nm NV height')
# plt.plot(xCrop, bzCropsNR[1,3][:, int(cropNum/2)], label='30 nm NV height')
# plt.plot(xCrop, bzCropsNR[1,4][:, int(cropNum/2)], label='40 nm NV height')
# plt.plot(xCrop, bzCropsNR[1,5][:, int(cropNum/2)], label='50 nm NV height')
# plt.legend(loc=3,borderaxespad=1,prop={'size':10})
# plt.title(u"Bz linecuts, right-handed Néel, 10 nm r0")

fig1, ax1 = plt.subplots()
plt.plot(xCrop, bzCropsNR[1,0][:, int(cropNum/2)],  label=u'right-handed Néel')
plt.plot(xCrop, bzCropsB[1,0][:, int(cropNum/2)], label='Bloch')
plt.plot(xCrop, bzCropsNL[1,0][:, int(cropNum/2)], label=u'left-handed Néel')
plt.legend(loc=4,borderaxespad=1,prop={'size':12})
plt.title("Bz linecuts, 10 nm FWHM, 10 nm NV height")
plt.ylabel(r"$B_z$ (T)")
plt.xlabel("x (nm)")
plt.savefig(savepath+'dwtype_fhwm_10nm.pdf',  bbox_inches='tight')

#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, bzCropsB[0,3][:, int(cropNum/2)], label='NV height = 30nm')
# plt.plot(xCrop, bzCropsB[0,4][:, int(cropNum/2)], label='NV height = 40nm')
# plt.plot(xCrop, bzCropsB[0,5][:, int(cropNum/2)], label='NV height = 50nm')
# plt.legend(loc=4,borderaxespad=1,prop={'size':10})
# plt.xlabel("x (nm)")
# plt.ylabel(r"$B_z$ (T)")
# plt.title(r"$B_z$ linecuts, Bloch, $r_0$ = 5 nm")
# plt.savefig(savepath+'bloch_cut_heights_far5nm.pdf',  bbox_inches='tight')
#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, bzCropsB[1,3][:, int(cropNum/2)], label='NV height = 30nm')
# plt.plot(xCrop, bzCropsB[1,4][:, int(cropNum/2)], label='NV height = 40nm')
# plt.plot(xCrop, bzCropsB[1,5][:, int(cropNum/2)], label='NV height = 50nm')
# plt.legend(loc=4,borderaxespad=1,prop={'size':10})
# plt.xlabel("x (nm)")
# plt.ylabel(r"$B_z$ (T)")
# plt.title(r"$B_z$ linecuts, Bloch, $r_0$ = 10 nm")
# plt.savefig(savepath+'bloch_cut_heights_far10nm.pdf',  bbox_inches='tight')
#
# fig1, ax1 = plt.subplots()
# plt.plot(xCrop, mz[cropRange[0]:cropRange[1], int(slen/2)])
# plt.xlabel("x (nm)")
# plt.ylabel(r"normalized $m_z$")
# plt.title(r"$m_z$ linecut, $r_0$ = 10 nm")
# plt.savefig(savepath+'mz_cut.pdf',  bbox_inches='tight')

fp.format_plots(plt, small=False)

plt.show()
