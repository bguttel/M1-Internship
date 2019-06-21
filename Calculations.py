# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:46:44 2019

@author: manip.batm
"""

import numpy as np

v= 2770
lambda_0 = 2e-6
f_c = 1.385*1e9
f = f_c - 0.1*1e9
k = (2*np.pi/v)*f
print('For f=%.3f [GHz], 2*pi/(k*lambda) = %.3f' %(f*1e-9,2*np.pi/(k*lambda_0)))