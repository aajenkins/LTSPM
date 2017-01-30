import numpy as np
import matplotlib.pyplot as plt
import format_plot as fp

# def analytic_m_calc(file,res):

pi = np.pi

filespec = "Msnotfixed"
file = "/Users/alec/UCSB/LTSPM/cofeb_analysis/ta/1760/stray_field_sim/radiusphi_"+filespec+".txt"
rphidata = np.loadtxt(file, delimiter=',')

rphiwrap0 = np.concatenate((np.add(rphidata[0,:-1],-2*pi),rphidata[0,:-1],np.add(rphidata[0,:-1],2*pi)))
rphiwrap1 = np.concatenate((rphidata[1,:-1],rphidata[1,:-1],rphidata[1,:-1]))

def phic(x, y):
	return np.arctan2(y, x) + pi

def rphi(phi):
	

def rphi_derivative(phi):


def dw_width(dw0, phi):
	1/np.cos(np.arctan( rphi(phi) / rphi ))


fig1, ax1 = plt.subplots()
plt.plot(rphi[0],rphi[1])
fp.format_plot(plt, 400, 400, 50, 50)

plt.show()