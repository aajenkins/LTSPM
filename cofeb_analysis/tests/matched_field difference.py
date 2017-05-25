# @Author: Jenkins Alec <alec>
# @Date:   2017-03-06T20:32:08-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-03-06T20:40:19-08:00



import matplotlib
import matplotlib.pyplot as plt
import numpy as np

path = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
filespec = 'Msfixed'
filespec = 'Msfixed'
cal_params_path = path+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

phi = cal_params['phi']
theta = cal_params['theta']
thetaError = cal_params['thetaError']

scannum = 1903

path = '/Users/alec/UCSB/cofeb_analysis_data/irmn/'
filespec = 'Msfixed'
cal_params_path = path+'cal_parameters_'+filespec+'.json'
with open(cal_params_path, 'r') as fread:
    cal_params = json.load(fread)

Ms = cal_params['Ms']
t = cal_params['t']
MstError = cal_params['MstError']
phi = cal_params['phi']
height = cal_params['height']
heightError = cal_params['heightError']
scanSize = cal_params['scanSize']

heights = [height - heighterr, height, height + heighterr]

basepath = "/Users/alec/UCSB/cofeb_analysis_data/irmn/stray_field_sim/"
savepath = "/Users/alec/UCSB/cofeb_analysis_data/irmn/stray_field_sim/"
