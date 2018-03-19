# @Author: Jenkins Alec <alec>
# @Date:   2018-03-18T14:16:17-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2018-03-19T14:32:01-07:00

import glob
import re
import numpy as np
import loadmat

class scan:

    def __init__(self, scanNum, xlen, ylen, basePath='/Users/alec/UCSB/scan_data/'):
        self.scanNum = scanNum
        self.dataFiles = self.get_data_files(basePath)
        self.param = loadmat.loadmat(self.dataFiles[0])['scan']['param']
        self.param['xlen'] = xlen
        self.param['ylen'] = ylen
        self.dataForward, self.dataReverse = self.load_data()


    def get_data_files(self, basePath):
        dirpath = basePath+str(self.scanNum)+'/'
        dataFiles = glob.glob(dirpath+'*.mat')
        numberedFiles = dict()
        for filePath in dataFiles:
            fileName = filePath.split('/')[-1]
            fileNum = int(re.search(r'[0-9]+.mat', fileName).group().split('.')[0])
            numberedFiles[fileNum] = fileName
        sortedDataFiles = [dirpath+numberedFiles[key] for key in sorted(numberedFiles.keys())]
        return sortedDataFiles

    def load_data(self):
        data = dict()
        dataForward = dict()
        dataReverse = dict()
        xlen = self.param['xlen']
        for key in loadmat.loadmat(self.dataFiles[0])['scan']['data'].keys():
            data[key] = []
        for file in self.dataFiles:
            fileData = loadmat.loadmat(file)['scan']['data']
            for key in data.keys():
                data[key].append(fileData[key])
        for key in data.keys():
            dataForward[key] = []
            dataReverse[key] = []
            for i in range(0, len(self.dataFiles), 2*xlen):
                dataForward[key].append(data[key][i : i + xlen])
                dataReverse[key].append(np.flip(data[key][i + xlen : i + 2*xlen], axis=0))
        return dataForward, dataReverse

    def set_NV_orientation(self, theta, phi):
        self.param['thetaNV'] = theta
        self.param['phiNV'] = phi
