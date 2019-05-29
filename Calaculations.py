# -*- coding: utf-8 -*-
"""
Created on Wed May 29 08:32:13 2019

@author: User
"""

#Constants

hbar = 6.5821195144*10e-16 #[eV*s]
kB = 8.61733035*10e-5 #[eV/K]
c = 299792458 #[m/s]

#Variables
T = 10*10e-3
omega = 2.3*10e9
v = 2770
WaveLength = v/omega

#Calculations

print((kB*T) / (hbar*omega))
