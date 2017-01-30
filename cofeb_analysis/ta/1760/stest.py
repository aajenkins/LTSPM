import numpy as np
import matplotlib.pyplot as plt
import format_plot as fp

path = '/Users/alec/UCSB/LTSPM/cofeb_analysis/ta/1760/stray_field_sim/stest.dat'

st = np.loadtxt(path)

fig1, ax1 = plt.subplots()
im1 = plt.imshow(st, cmap='gray', interpolation='nearest')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 500, 500, 0, 50)
plt.show()