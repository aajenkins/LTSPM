# @Author: Jenkins Alec <alec>
# @Date:   2017-08-01T18:26:53-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-01T18:27:12-07:00



import matplotlib.pyplot as plt
from . import format_plots_tkagg as fp

def aimshow(*args, **kwargs):
    plt.imshow(*args, **kwargs)
    fp.format_plot(plt, fwidth=600, fheight=600)
