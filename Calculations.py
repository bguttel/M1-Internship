# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:46:44 2019

@author: manip.batm
"""



import numpy as np
import matplotlib.pyplot as plt
import time
start_time = time.time()

def abs2(x):
    return x.real**2 + x.imag**2

def Z2Y_matrix(Z):
    Y = np.zeros((2,2), dtype=np.complex) #Y-matrix
    det_Z = Z[0,0]*Z[1,1] - Z[0,1]*Z[1,0]
    Y[0,0] = Z[1,1]/det_Z
    Y[0,1] = -Z[0,1]/det_Z
    Y[1,0] = -Z[1,0]/det_Z
    Y[1,1] = Z[0,0]/det_Z
    return Y

def Y2S_matrix(Y,Y0):
    """S-matrix of a delay-line using Y-matrix"""
###  Based on Morgan's book p.403  
#    S[0,0] = ((1-R1*Y[0,0])*(1+R2*Y[1, 1])+R1*R2*np.square(Y[0,1]))/((1+R1*Y[0,0])*(1+R2*Y[1,1])-R1*R2*np.square(Y[0,1]))
#    S[1,1] = ((1-R2*Y[1,1])*(1+R1*Y[0,0])+R2*R1*np.square(Y[1,0]))/((1+R2*Y[1,1])*(1+R1*Y[0,0])-R1*R2*np.square(Y[1,0]))
#    S[0,1] = S[1,0] = 2*Y[0,1]*np.sqrt(R1*R2)/(R1*R2*np.square(Y[0,1])-(1+Y[0,0]*R1)*(1+Y[1,1]*R2))

###  Based on Pozar's book p.192
#    S = np.zeros((2,2), dtype=np.complex) #P-matrix
#    Delta_Y = (Y[0,0]+Y0)*(Y[1,1]+Y0)-Y[0,1]*Y[1,0]
#    S[0,0] = ((Y0-Y[0,0])*(Y0+Y[1,1])+Y[0,1]*Y[1,0])/Delta_Y
#    S[1,1] = ((Y0+Y[0,0])*(Y0-Y[1,1])+Y[0,1]*Y[1,0])/Delta_Y
#    S[0,1] = S[1,0] = (-2*Y[0,1]*Y0)/Delta_Y
    
    S = np.dot(np.linalg.inv(np.linalg.inv(Y)+(1/Y0)*np.identity(2)), np.linalg.inv(Y)-(1/Y0)*np.identity(2))
    print(abs2(S[0,1])+abs2(S[0,0])) 
    return S

print(2770/0.42)
                        
                        

#epsilon_0 = 8.8541878e-12 #[F/m]
#j = np.complex(0,1)
#v_f= 3845
#v= v_f
#lambda_0 = 800e-9
#f_0 = 4.8e9 #design frequency
#v_0 = lambda_0 * f_0
#N = 135
#f_c = v/lambda_0
#f = f_c - 0.1*1e9
#k = (2*np.pi/v)*f
##print('For f=%.3f [GHz], 2*pi/(k*lambda) = %.3f' %(f*1e-9,2*np.pi/(k*lambda_0)))
####Far from center frequency S-matrix calculation (simple circuit approximation)
#
#num_frequencies = 10
#frequencies_width = 2* (7 * 1e9 - f_c) #Desired total width of frequencies simulated, <=f_c
#frequencies = np.linspace(f_c-(frequencies_width/2),f_c+(frequencies_width/2),num_frequencies) #sent frequencies
#
#S_far = np.zeros((num_frequencies,2,2), dtype=np.complex)
#Y_far = np.zeros((num_frequencies,2,2), dtype=np.complex)
#Z_far = np.zeros((num_frequencies,2,2), dtype=np.complex)
#
#for i,f in enumerate(frequencies):
#    omega = f*2*np.pi
#    R_x = 15
#    C_s = 0.9e-12
#    Y0 = 1/50
#        
#    Z_far[i,0,0] = Z_far[i,1,1] = 1/Y0
#    Z_far[i,0,1] = Z_far[i,1,0] = 1/(j*omega*C_s) +2*R_x
#    Y_far[i,:,:] = Z2Y_matrix(Z_far[i,:,:])
#
##    Y0 = np.sqrt((np.real(Y_far[i,0,0])+j*omega*C_s)/(2*R_x))
#    S_far[i,:,:] = Y2S_matrix(Y_far[i,:,:], Y0)
#    
##plot of behavior far from center frequency (capacitive stray effect)
#fig3, axs_far = plt.subplots(2)
#axs_far[0].plot(frequencies*1e-9,20*np.log10(abs2(S_far[:,1,0])), '-')
#axs_far[1].plot(frequencies*1e-9,20*np.log10(abs2(S_far[:,0,0])), '-')
#axs_far[1].set_ylabel(r'$|S_{11}|^2$ [dB]')
#axs_far[0].set_ylabel(r'$|S_{21}|^2$ [dB]')
#axs_far[1].set_xlabel(r'Frequency [GHz]')
#axs_far[0].set_xlabel(r'Frequency [GHz]')
#fig3.suptitle(r'$S_{far}$' '\n' r'$Z_0=%.1f$' %(1/Y0))
#
#
#plt.show()
print("--- %s seconds ---" % (time.time() - start_time))