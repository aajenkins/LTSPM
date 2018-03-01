# @Author: Jenkins Alec <alec>
# @Date:   2017-07-09T16:46:45-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-30T12:52:28-07:00


import numpy as np
import matplotlib.pyplot as plt
import json

import plotting.format_plots_tkagg as fp
from PIL import Image
import calc_NV_field as cNVf
import load_scan as lscan
import linecut

pi = np.pi

def plot_stray_field(scannum, helicities = [0, 90, 180]):

    path = '/Users/alec/UCSB/cofeb_analysis_data/ta/'
    scan_params_path = path+str(scannum)+'/'+'scan_parameters.json'
    field_path = path+str(scannum)+'/stray_field_sim/'
    datapath = path+'1760/'

    with open(scan_params_path, 'r') as fread:
        scan_params = json.load(fread)

    phi = scan_params['phi']
    theta = scan_params['theta']
    xres = scan_params['xres']
    xcenter = scan_params['xcenter']
    ycenter = scan_params['ycenter']
    scanSize = (1e6)*scan_params['scanSize'] # convert to microns

    errnames = ["lower", "mean", "upper"]

    ffdata = lscan.load_ff('/Users/alec/UCSB/scan_data/'+str(scannum)+'-esrdata/fitdata.txt',xres,xres,maxfgrad=20)

    scd = np.empty((len(helicities)), dtype=object)
    vcd = np.empty_like(scd)
    meff = np.empty_like(scd)
    bx = np.empty((len(helicities), len(errnames)), dtype=object)
    by = np.empty_like(bx)
    bz = np.empty_like(bx)
    bNV = np.empty_like(bx)
    fieldCutNV = np.empty_like(bx)

    for j in range(len(helicities)):
        scd[j] = np.loadtxt(field_path+'scd_'+str(helicities[j])+str(scannum)+'.txt', delimiter=',')
        vcd[j] = np.loadtxt(field_path+'vcd_'+str(helicities[j])+str(scannum)+'.txt', delimiter=',')
        meff[j] = np.loadtxt(field_path+'meff_'+str(helicities[j])+str(scannum)+'.txt', delimiter=',')
        for i in range(len(errnames)):
            bx[j,i] = (1e4)*np.loadtxt(field_path+'h'+str(helicities[j])+'_x_'+errnames[i]+'_'+str(scannum)+'.txt', delimiter=',')
            by[j,i] = (1e4)*np.loadtxt(field_path+'h'+str(helicities[j])+'_y_'+errnames[i]+'_'+str(scannum)+'.txt', delimiter=',')
            bz[j,i] = (1e4)*np.loadtxt(field_path+'h'+str(helicities[j])+'_z_'+errnames[i]+'_'+str(scannum)+'.txt', delimiter=',')
            bNV[j,i] = cNVf.calc_NV_field(bx[j,i], by[j,i], bz[j,i], theta, phi)

    slen = len(bx[0,0])

    cutSize = 2.2
    phinum = 8
    philist = np.linspace(0, pi, phinum, endpoint=False)

    ffdataCut = np.empty((phinum), dtype=object)
    ffdataCutErr = np.empty((phinum), dtype=object)
    fieldCutNV = np.empty((phinum, len(helicities), len(errnames)), dtype=object)

    for k in range(phinum):
        ffdataCut[k] = linecut.linecut(ffdata[0], scanSize, cutSize, philist[k], xcenter, ycenter)
        ffdataCutErr[k] = linecut.linecut(ffdata[1], scanSize, cutSize, philist[k], xcenter, ycenter)
        for j in range(len(helicities)):
            for i in range(len(errnames)):
                fieldCutNV[k,j,i] = linecut.linecut(bNV[j,i], scanSize, cutSize, philist[k])


    #---------------- PLOTS ------------------------------------------
    #-----------------------------------------------------------------

    savepath = '/Users/alec/UCSB/papers/tacofeb/figures/'

    plt.close('all')

    if(len(helicities)==4):
        fig, ax = plt.subplots(figsize=(4,4))
        im = plt.imshow(bNV[3,1], interpolation='nearest', cmap='bone')
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        plt.subplots_adjust(left=0.0, bottom=0.0, right=0.9, top=1.0, wspace=0, hspace=0)
        plt.savefig(savepath+'BNV_bestFitHelicity_'+str(scannum)+'.pdf', format='pdf')

    fig, ax = plt.subplots(figsize=(4,4))
    im = plt.imshow(ffdata[0], interpolation='nearest', cmap='bone')
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.subplots_adjust(left=0.0, bottom=0.0, right=0.9, top=1.0, wspace=0, hspace=0)
    plt.savefig(savepath+'BNV_data_'+str(scannum)+'.pdf', format='pdf')

    fig, axes = plt.subplots(ncols=2, nrows=int(phinum/2), sharex=True, sharey=True, figsize=(6,7.5))
    for j in range(0,2):
        for i in range(0,int(phinum/2)):
            axes[i,j].errorbar(ffdataCut[int(i+(phinum/2)*j)][0], ffdataCut[int(i+(phinum/2)*j)][1], color='#000000', fmt='.', linewidth=1.0, label='data')
            axes[i,j].plot(fieldCutNV[0,0,1][0], fieldCutNV[int(i+(phinum/2)*j),0,1][1], linewidth=2.0, label=u'right-handed Néel')
            axes[i,j].plot(fieldCutNV[0,1,1][0], fieldCutNV[int(i+(phinum/2)*j),1,1][1], linewidth=2.0, label="Bloch")
            axes[i,j].plot(fieldCutNV[0,2,1][0], fieldCutNV[int(i+(phinum/2)*j),2,1][1], linewidth=2.0, label=u'left-handed Néel')
            if(len(fieldCutNV[0,:,1])==4):
                axes[i,j].plot(fieldCutNV[0,3,1][0], fieldCutNV[int(i+(phinum/2)*j),3,1][1], linewidth=2.0, label=r'$\psi_h$ = '+str(helicities[3])+u'°')
            axes[i,j].get_yaxis().set_visible(False)
            axes[i,j].get_xaxis().set_visible(False)
            axes[i,j].text(0.04,0.86,r'$\phi$ = '+'{:d}{}{:d}'.format(int(i+(phinum/2)*j),r'$\pi$/',phinum),
                horizontalalignment='left', verticalalignment='center',
                transform=axes[i,j].transAxes, fontsize=10)
    axes[int((phinum/2)-1),1].get_yaxis().set_visible(True)
    axes[int((phinum/2)-1),1].get_xaxis().set_visible(True)
    axes[int((phinum/2)-1),1].yaxis.tick_right()
    axes[int((phinum/2)-1),1].yaxis.set_label_position("right")
    axes[int((phinum/2)-1),1].set_xlabel(r'r ($\mu$m)')
    axes[int((phinum/2)-1),1].set_ylabel(r'B$\mathrm{_{NV}}$ (G)')
    axes[int((phinum/2)-1),0].legend(bbox_to_anchor=(0.0, -0.74), loc=3, borderaxespad=0., frameon=False, prop={'size':10})
    plt.ylim([0,32])
    plt.subplots_adjust(left=0.0, bottom=0.15, right=0.88, top=1.0, wspace=0, hspace=0)
    plt.savefig(savepath+'BNV_linecuts_'+str(scannum)+'.pdf', format='pdf')


    # fp.format_plots(plt, small=False, tight=False)

    plt.show()

if __name__ == "__main__":
    import sys
    if (len(sys.argv) == 2):
        plot_stray_field(int(sys.argv[1]))
    elif (len(sys.argv) == 3):
        plot_stray_field(int(sys.argv[1]), helicities=np.array(eval(sys.argv[2])))
    else:
        print('enter scan number')
