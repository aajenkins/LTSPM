# @Author: Jenkins Alec <alec>
# @Date:   2017-07-09T13:26:11-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-07-09T13:31:38-07:00



from ctypes import cdll
from ctypes import c_char_p

lib = cdll.LoadLibrary('test.so')
hello = lib.fun
hello()
