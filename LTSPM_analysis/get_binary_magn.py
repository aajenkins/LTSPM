# @Author: Jenkins Alec <alec>
# @Date:   2017-02-03T18:34:58-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-04T20:53:54-08:00



import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import imread

def get_binary_magn(impath):
    domains = imread(impath, flatten=True)
    domains = np.add(np.multiply(2/255,domains),-1)

    return domains
