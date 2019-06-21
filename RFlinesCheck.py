# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 08:49:58 2019

@author: manip.batm
"""

import numpy as np
#import h5py
import matplotlib.pyplot as plt

frequency = []
S11_0 = []
S12_0 = []
S21_0 = []
S22_0 = []

###Get calibration data
with open('C:/Users/manip.batm/Desktop/DClines/RF/17-06-19/0.s2p', "r") as f:
    lines = f.readlines()
     
    # Loop through all lines, ignoring header.
    # Add last element to list (i.e. the process name)
    for l in lines[8:]:
        frequency.append(float(l.split()[0]))
        S11_0.append(float(l.split()[1]))
        S21_0.append(float(l.split()[3]))
        S12_0.append(float(l.split()[5]))
        S22_0.append(float(l.split()[7]))
frequency = np.asarray(frequency)
S11_0 = np.asarray(S11_0)
S12_0 = np.asarray(S12_0)
S21_0 = np.asarray(S21_0)
S22_0 = np.asarray(S22_0)

###Actual data            
RFs = ['1','2','3','4']
colors = ['b','k','r','g']
S21_1 = []
S21_2 = []
S21_3 = []
S21_4 = []
S21s = [S21_1,S21_2,S21_3,S21_4]


for i,line in enumerate(RFs):
    S11 = []
    S12 = []
    S21 = []
    S22 = []

# =============================================================================
# ###Plot all graphs for every RF line
#     
    with open('C:/Users/manip.batm/Desktop/DClines/RF/17-06-19/'+line+'.s2p', "r") as f:
         lines = f.readlines()
      
         # Loop through all lines, ignoring header.
         # Add last element to list (i.e. the process name)
         for l in lines[8:]:
             S11.append(float(l.split()[1]))
             S21.append(float(l.split()[3]))
             S12.append(float(l.split()[5]))
             S22.append(float(l.split()[7])) 
#             
    S11 = np.asarray(S11)
    S12 = np.asarray(S12)
    S21 = np.asarray(S21)
    S22 = np.asarray(S22)
    S21s[i] = S21-S21_0
#     
#    figure, ax = plt.subplots()
#    ax.plot(frequency*1e-9,S11+S11_0, label='S11')
#    ax.plot(frequency*1e-9,S12+S12_0, label='S12')
#    ax.plot(frequency*1e-9,S21+S21_0, label='S21')
#    ax.plot(frequency*1e-9,S22+S22_0, label='S22')
#    ax.set_xscale('log')
#    ax.set_title('RF line no.'+line)
#    ax.set_xlabel('Frequency [GHz]')
#    ax.set_ylabel('Attenuation [dB]')
#    ax.legend()   
# =============================================================================
###Plot S21 for all RF lines
    
    plt.plot(frequency*1e-9,S21s[i], colors[i], label='RF'+RFs[i])
    
plt.title('S21')
plt.xlabel('Frequency [GHz]')
plt.xscale('log')
plt.yticks(np.arange(-50,0,10))
plt.ylabel('Attenuation [dB]')
plt.ylim(-20,0)
plt.yticks(np.arange(-20, 2.5, step=2.5))
plt.legend()
# =============================================================================

plt.show()