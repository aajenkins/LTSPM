# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 21:10:29 2016

@author: alec
"""

import numpy as np

info = open('/Users/alec/UCSB/scan_data/809/000809.info','r')

lines = info.readlines()

size = float(lines[1][6])*5

print(size)
