# @Author: Jenkins Alec <alec>
# @Date:   2017-08-01T17:31:05-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-08-01T18:35:56-07:00



import matplotlib.pyplot as plt
from . import format_plots_tkagg as fp

def aplot(*args, **kwargs):
    fig, ax = plt.subplots()
    plt.plot(*args, **kwargs)
    fp.format_plot(plt)

def aimshow(*args, **kwargs):
    fig, ax = plt.subplots()
    plt.imshow(*args, **kwargs)
    fp.format_plot(plt, fwidth=600, fheight=600)
