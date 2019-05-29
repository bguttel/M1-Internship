# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:11:15 2019

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt

#Gives the equivalent resistance of two parallel resistors
def ParaRes(R1,R2):
    R3 = 1/( (1/R1) + (1/R2) )
    return R3

#stretch an array to a desire length by adding linear interpolations between values
def Stretch(a,L):
    b = []
#    d - length of regular blocks, r - reminder, number of longer blocks
    d = int((L-1)/(len(a)-1))
    r = (L-1)-((len(a)-1)*d)
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

a=np.array([1,2,3,4,4,10,15])   
L = 17
b = Stretch(a,L)
#print(a)
#print (b)
plt.plot(np.linspace(0,L,len(a)),a, 'ro')
plt.plot(np.linspace(0,L,L),b, 'b-')
print(len(b)-L)