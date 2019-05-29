# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:36:30 2019

@author: User
"""


import numpy as np
import matplotlib.pyplot as plt

#Gives the equivalent resistance of two parallel resistors
def ParaRes(R1,R2):
    R3 = 1/( (1/R1) + (1/R2) )
    return R3
 
#print(ParaRes(5,5))

#Wires Resistances
R1 = 170.5
R15 = 165.6
R16 = 231.3
R17 = 233.1
R21 = 232.2
R22 = 232.1
R23 = 166.5
R24 = 242.3
R26 = 167.2
R40 = 166.3
R42 = 166.5

#Disconnected Wire
dis = float(10**20)

RA = ParaRes(R24,R1) + R26
RB = ParaRes(R22,R42) + R23
print('In this case the total resistance is' ,round(RA,2))
