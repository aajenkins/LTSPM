# @Author: Jenkins Alec <alec>
# @Date:   2017-07-27T21:26:18-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-28T09:19:24-07:00

import numpy as np
from scipy import ndimage


def linecut(data, size, cutSize, angle, xcenter=0, ycenter=0):

    dlen = len(data)
    dres = size/dlen
    cutlen = np.round(cutSize/(2*dres))

    if (xcenter == 0 and ycenter == 0):
        xcenter = int(dlen/2)
        ycenter = int(dlen/2)

    s = np.arange(-dres*cutlen, dres*(cutlen-(1/2)), dres)

    cutx1, cuty1 = xcenter - cutlen*np.cos(angle), ycenter - cutlen*np.sin(angle)
    cutx2, cuty2 = xcenter + cutlen*np.cos(angle), ycenter + cutlen*np.sin(angle)
    x, y = np.linspace(cutx1, cutx2, 2*cutlen, endpoint=False), np.linspace(cuty1, cuty2, 2*cutlen, endpoint=False)

    dataCut = ndimage.map_coordinates(np.transpose(data), np.vstack((x,y)), order=1)

    return [s, dataCut]
