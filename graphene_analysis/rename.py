# @Author: Jenkins Alec <alec>
# @Date:   2018-03-12T19:01:02-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2018-03-18T11:02:03-07:00



import os
import glob
import re

dirpath = '/Users/alec/UCSB/scan_data/2074/'
filepaths = glob.glob(dirpath+'scanData*.mat')
filenames = [f.split('/')[-1] for f in filepaths]

for i in range(len(filenames)):
    # fnum = re.search(r'scanData__[0-9]+', filenames[i]).group()
    # fnum = int(re.sub(r'scanData__0+', '', fnum))
    # newfnum = str(int(fnum-243))
    # newfilename = 'scanData2070_'+newfnum+'.mat'
    newfilename = re.sub(r'2070', '2074', filenames[i])
    os.rename(dirpath+filenames[i], dirpath+newfilename)
