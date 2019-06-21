# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:11:15 2019

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt

#Gives the equivalent resistance of two parallel resistors
def ParaRes(*Rs):
    R_rec = 0
    for R in Rs:
        R_rec = R_rec+(1/R)
    R_eqv = 1/R_rec
    return R_eqv

#stretch an array to a desire length by adding linear interpolations between values
def Stretch(a,L):
    b = []
#    d - length of regular blocks, r - reminder, number of longer blocks
    d = int(L/(len(a)-1)-1)
    r = L-1-((len(a)-1)*d)
    for i in range (len(a)-1):
        if r==0:
            c = np.interp(np.linspace(0,d-1,d) , [0,d], [a[i],a[i+1]])
            b = np.concatenate((b, c))
        else :
            c = np.interp(np.linspace(0,d,d+1) , [0,d+1], [a[i],a[i+1]])
            b = np.concatenate((b, c))
            r = r-1
    b = np.append(b,a[-1])
    return b