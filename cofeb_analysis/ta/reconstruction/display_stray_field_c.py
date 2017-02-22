# @Author: Jenkins Alec <alec>
# @Date:   2017-02-09T16:32:53-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-15T20:18:26-08:00



import numpy as np
import matplotlib.pyplot as plt
import format_plot as fp

pi = np.pi

path = '/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/'
hz = np.loadtxt(path+'stray_field_test_c.txt', delimiter=',')

plt.close('all')

fig1, ax1 = plt.subplots()
im1 = plt.imshow(hz, cmap='jet')
fig1.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
fp.format_plot(plt, 400, 400, 50, 50)

plt.show()
