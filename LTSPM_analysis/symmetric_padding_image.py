# @Author: Jenkins Alec <alec>
# @Date:   2017-02-03T15:13:39-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-03T15:16:17-08:00

# returns original array with mirror padding of original array size

import numpy as np

def symmetric_padding_image(image):
    image1 = np.flipud(np.fliplr(image))
    image2 = np.flipud(image)
    image3 = np.flipud(np.fliplr(image))
    image4 = np.fliplr(image)
    image5 = np.fliplr(image)
    image6 = np.flipud(np.fliplr(image))
    image7 = np.flipud(image)
    image8 = np.flipud(np.fliplr(image))

    imagetop = np.concatenate((image1,image2,image3), axis=1)
    imagemiddle = np.concatenate((image4,image,image5), axis=1)
    imagebottom = np.concatenate((image6,image7,image8), axis=1)

    return np.concatenate((imagetop, imagemiddle, imagebottom),
                                    axis=0)
