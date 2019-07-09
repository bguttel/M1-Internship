# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:11:06 2019

@author: manip.batm
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
import time
start_time = time.time()

def legendre_function(s,x):
    a = 1
    Ps_x = 1
    for m in range(1,50):
        a = a*(m-1-s)*(m+s)*(1-x)/(2*m**2)
        Ps_x += a
    return Ps_x 

def rho_f(k,a,p,epsilon_infty):
    """Analytical form of the Fourier transform of the elemtal charge density"""
    m = np.floor((k*p)/(2*np.pi))
    s = ((k*p)/(2*np.pi))-m
    Delta = np.pi*a/p
    return epsilon_infty * 2 *np.sin(np.pi*s)* (legendre_function(m,np.cos(Delta))/legendre_function(-s,-np.cos(Delta)))

def rho_e(N,S_e,k,a,p,epsilon_infty):
    """Analytical form of the Fourier transform of a uniform tranducer's charge density. Based on Morgan's book, p.141"""
    if S_e == 2 or S_e == 3: #number of electrodes per cell
        Abar_1 = 1
    elif S_e == 4:
        Abar_1 = 2*np.cos(k*p/2)
    Abar_N = np.sin(N*k*p*S_e/2)/np.sin(k*p*S_e/2)
    return Abar_1 * Abar_N *rho_f(k,a,p,epsilon_infty)

def FWHM(x, y):
    """Gives the full width of half max of a nice-enough function"""
    half = (min(y)+max(y))/2
    i=0
    while y[i]<half:
        i += 1
    left_idx = i
    i=1
    while y[-i]<half:
        i += 1
    right_idx = i
    return x[-right_idx]-x[left_idx]

def DART_to_iDART(P):
    """Invert ports 1 and 2"""
#    P_inverted = np.zeros((3,3), dtype=np.complex) #P-matrix
    P_inverted = np.copy(P)
    P_inverted[0,0],P_inverted[1,1] = P[1,1],P[0,0]
    P_inverted[0,1],P_inverted[1,0] = P[1,0],P[0,1]
    P_inverted[0,2],P_inverted[1,2] = P[1,2],P[0,2]
    P_inverted[2,0],P_inverted[2,1] = P[2,1],P[2,0]
    
    return P_inverted    

def P2Y_matrix(Pa,Pb):
    """Y-matrix of a delay-line using P-matrices"""
    Y = np.zeros((2,2), dtype=np.complex) #Y-matrix
    Y[0,0] = Pa[2,2] - 2*Pb[0,0]*np.square(Pa[0,2])/(1-Pa[0,0]*Pb[0,0])
    Y[1,1] = Pb[2,2] - 2*Pa[0,0]*np.square(Pb[0,2])/(1-Pa[0,0]*Pb[0,0])
    Y[0,1] = Y[1,0] = -2*Pa[0,2]*Pb[0,2]/(1-Pa[0,0]*Pb[0,0])
       
    return Y

def Z2Y_matrix(Z):
    Y = np.zeros((2,2), dtype=np.complex) #Y-matrix
    det_Z = Z[0,0]*Z[1,1] - Z[0,1]*Z[1,0]
    Y[0,0] = Z[1,1]/det_Z
    Y[0,1] = -Z[0,1]/det_Z
    Y[1,0] = -Z[1,0]/det_Z
    Y[1,1] = Z[0,0]/det_Z
    return Y

def Y2S_matrix(Y,R1,R2):
    """S-matrix of a delay-line using Y-matrix"""
###  Based on Morgan's book p.403  
#    S[0,0] = ((1-R1*Y[0,0])*(1+R2*Y[1, 1])+R1*R2*np.square(Y[0,1]))/((1+R1*Y[0,0])*(1+R2*Y[1,1])-R1*R2*np.square(Y[0,1]))
#    S[1,1] = ((1-R2*Y[1,1])*(1+R1*Y[0,0])+R2*R1*np.square(Y[1,0]))/((1+R2*Y[1,1])*(1+R1*Y[0,0])-R1*R2*np.square(Y[1,0]))
#    S[0,1] = S[1,0] = 2*Y[0,1]*np.sqrt(R1*R2)/(R1*R2*np.square(Y[0,1])-(1+Y[0,0]*R1)*(1+Y[1,1]*R2))

###  Based on Pozar's book p.192
    S = np.zeros((2,2), dtype=np.complex) #P-matrix
    Y0 = np.sqrt(R1*R2)
    Delta_Y = (Y[0,0]+Y0)*(Y[1,1]+Y0)-Y[0,1]*Y[1,0]
    S[0,0] = ((Y0-Y[0,0])*(Y0+Y[1,1])+Y[0,1]*Y[1,0])/Delta_Y
    S[1,1] = ((Y0+Y[0,0])*(Y0-Y[1,1])+Y[0,1]*Y[1,0])/Delta_Y
    S[0,1] = S[1,0] = -2*Y[0,1]*Y0/Delta_Y
        
    return S

def S21(Pa,Pb,Y0):
    """Transmission coefficient of two P-matrices"""
    Delta_Y = (2*np.square(Pa[0,2])*(2*np.square(Pb[0,2])+Pb[0,0]*(Pb[2,2]+Y0))+
               (Pa[2,2]+Y0)*(2*Pa[0,0]*np.square(Pb[0,2])-Pb[2,2]-Y0+Pa[0,0]*Pb[0,0]*(Pb[2,2]+Y0)))
    return -4*Pa[0,2]*Pb[0,2]*Y0/Delta_Y
    
def abs2(x):
    return x.real**2 + x.imag**2
    
    
### Natural constants
epsilon_0 = 8.8541878e-12 #[F/m]
j = np.complex(0,1)

###Material (GaAs) constants
epsilon_infty = 56*epsilon_0 #effective permittivity
velocities_difference_percent = 0.019 #(v_f-v_m)/v_f
Gamma_s = velocities_difference_percent/epsilon_infty #coupling constant
v_f = 3865 #[m/s] free velocity in the material
v_m = (1-velocities_difference_percent)*v_f #[m/s] velocity on metal
R = np.complex(0,0.04) #reflection constant
T = np.sqrt(1-abs2(R)) #transmission constant
C_S = 3.575e-09 #[F/m] Effective capacitance 
K_squared = velocities_difference_percent*2 #Piezo-electric coupling coefficient

#lambda_0s = np.asarray([0.5,0.75,1,2])*1e-6 #wavelength
lambda_0 = 800*1e-9 #wavelength
#Ns = np.array([90,110,130]) #total number of periods
N = 135 #total number of periods


#Half_widths = np.zeros((lambda_0s.size,Ns.size))
#Max_amplitudes = np.zeros((lambda_0s.size,Ns.size))

###Design constants in S.I units 
f_c = v_f/lambda_0 #center frequency 
W = 21*1e-6 #aperture
L = N*lambda_0 #total langth
p = lambda_0 /4 #pitch
eta = 0.5 #mettalization ratio
a = p*eta #mettalized electrode width
omega_c = 2*np.pi*f_c #center frequency
k_c = omega_c / v_f



num_frequencies = 300
frequencies_width = 1 * 1e9 #Desired total width of frequencies simulated, <=f_c
frequencies = np.linspace(f_c-(frequencies_width/2),f_c+(frequencies_width/2),num_frequencies) #sent frequencies



"""P-matrix calculations using RAM and COM analyses"""
P_DART_RAM = np.zeros((num_frequencies,3,3), dtype=np.complex) #RAM DART P-matrix
P_DART_COM = np.zeros((num_frequencies,3,3), dtype=np.complex) #COM DART P-matrix
P_iDART_RAM = np.zeros((num_frequencies,3,3), dtype=np.complex) #RAM inverted DART P-matrix
P_iDART_COM = np.zeros((num_frequencies,3,3 ), dtype=np.complex) #COM inverted DART P-matrix
P_double = np.zeros((num_frequencies,3,3), dtype=np.complex) #double IDT P-matrix
S_far = np.zeros((num_frequencies,2,2), dtype=np.complex)
Y_far = np.zeros((num_frequencies,2,2), dtype=np.complex)
Z_far = np.zeros((num_frequencies,2,2), dtype=np.complex)
c = np.zeros(N+1, dtype=np.complex) #forward wave
b = np.zeros(N+1, dtype=np.complex) #backwards wave
I = np.zeros(N+1, dtype=np.complex) #current

Ga_DART_RAM = np.zeros(num_frequencies)
Ga_DART_COM = np.zeros(num_frequencies)
Ga_double = np.zeros(num_frequencies)
#        c12 = np.zeros(num_frequencies, dtype=np.complex)
#Dram = np.zeros(num_frequencies)
Dcom = np.zeros(num_frequencies)
Ga = np.zeros(num_frequencies)
conservation_check = np.zeros(num_frequencies, dtype=np.complex)

###     RAM analysis
for i,f in enumerate(frequencies):
    omega = f*2*np.pi
        
    Delta_v_e =  ((v_m-v_f)/2)*(1+legendre_function((omega*p)/(2*np.pi*v_f),-np.cos(np.pi*a/p))/
              legendre_function(((omega*p)/(2*np.pi*v_f))-1,-np.cos(np.pi*a/p))) #[m/s] change in velocity due to electrical loading
    v = v_f + Delta_v_e #[m/s] effective wave velocity in the material, Morgan's book p.237


    k = omega/v #wavenumber
    
    c[0] = 0
    b[0] = 1
    for n in range(1,N+1):
        c[n] = (1/T)*c[n-1]*np.exp(-j*k*lambda_0) + (R/T)*b[n-1]
        b[n] = -(R/T)*c[n-1] + (1/T)*b[n-1]*np.exp(j*k*lambda_0)
        I[n] = -j*omega*W*rho_f(k,a,p,epsilon_infty)*(c[n]*np.exp(-j*k*(p/2))+b[n]*np.exp(j*k*(p/2)))*np.sqrt((2*Gamma_s)/(omega*W))
    I_total = np.sum(I)
        
    P_DART_RAM[i,1,1] = c[N]/b[N]
    P_DART_RAM[i,0,1] = b[0]/b[N]
    P_DART_RAM[i,2,1] = I_total/b[N]
    P_DART_RAM[i,1,2] = -P_DART_RAM[i,2,1]/2
    
    c[N] = 1
    b[N] = 0
    for n in range(1,N+1):
        c[-(n+1)] = (1/T)*c[-n]*np.exp(j*k*lambda_0) - (R/T)*b[-n]
        b[-(n+1)] = (R/T)*c[-n] + (1/T)*b[-n]*np.exp(-j*k*lambda_0)
        I[-(n+1)] = -j*omega*W*rho_f(k,a,p,epsilon_infty)*(c[-(n+1)]*np.exp(-j*k*(p/2))+b[-(n+1)]*np.exp(j*k*(p/2)))*np.sqrt((2*Gamma_s)/(omega*W))
    I_total = np.sum(I)
        
    P_DART_RAM[i,0,0] = b[0]/c[0]
    P_DART_RAM[i,1,0] = c[N]/c[0]
    P_DART_RAM[i,2,0] = I_total/c[0]
    P_DART_RAM[i,0,2] = -P_DART_RAM[i,2,0]/2

#            Dram[i] = np.abs(P[0,2]/P[1,2])
    Ga_DART_RAM[i] = np.abs(P_DART_RAM[i,0,2])**2 + np.abs(P_DART_RAM[i,1,2])**2
#           conservation_check[i] = (P[0,0]*np.conj(P[0,2])+P[0,1]*np.conj(P[1,2]+P[0,2]))
 
###         COM analysis
    c12 = np.conj( -(R/lambda_0)*np.exp(-2*j*k*(lambda_0*(3/4))) )
    aT = rho_f(k,a,p,epsilon_infty)*np.sqrt(omega*W*Gamma_s/2)
    alpha1 = np.conj(-j*(aT/lambda_0)*np.exp(-j*k*lambda_0*(3/8)))
    theta = np.angle(alpha1)
    phi = np.angle(c12)
    P_DART_COM[i,0,0] = -(np.conj(c12)/np.abs(c12))*np.tanh(np.abs(c12)*L)
    P_DART_COM[i,1,1] = -(c12/np.abs(c12))*np.tanh(np.abs(c12)*L)*np.exp(-2*j*k_c*L)
    P_DART_COM[i,0,1] = P_DART_COM[i,1,0] = np.exp(-j*k*L)/np.cosh(np.abs(c12)*L)
    delta = k - k_c
    s = np.sqrt(np.complex(delta**2 - abs2(c12)))
    K1 = (np.conj(alpha1)*c12-j*delta*alpha1)/(s**2)
    K2 = (np.conj(c12)*alpha1+j*delta*np.conj(alpha1))/(s**2)
    D = s*np.cos(s*L) + j*delta*np.sin(s*L)
    P_DART_COM[i ,2,0] = (2*np.conj(alpha1)*np.sin(s*L) - 2*s*K2*(np.cos(s*L)-1))/D
    P_DART_COM[i,2,1] = np.exp(-j*k_c*L)*((-2*alpha1*np.sin(s*L) - 2*s*K1*(np.cos(s*L)-1))/D)
    P_DART_COM[i,0,2] = P_DART_COM[i,2,0]/(-2)
    P_DART_COM[i,1,2] = P_DART_COM[i,2,1]/(-2)
    Ga_DART_COM[i] = abs2(P_DART_COM[i,0,2]) + abs2(P_DART_COM[i,1,2])
#    Dcom[i] = (np.abs(((1/np.tanh(np.abs(c12)*L/2))+np.exp(j*(2*theta-phi)))/((1/np.tanh(np.abs(c12)*L/2))-np.exp(-j*(2*theta-phi)))))**N

###Far from center frequency S-matrix calculation (simple circuit approximation)
    R_x = 15
    C_s = 0.9e-12
    Z_far[i,0,0] = Z_far[i,1,1] = R_x
    Z_far[i,0,1] = Z_far[i,1,0] = 1/(j*omega*C_s) +2*R_x
    Y_far[i,:,:] = Z2Y_matrix(Z_far[i,:,:])    
    S_far[i,:,:] = Y2S_matrix(Y_far[i,:,:], np.sqrt((np.real(Y_far[i,0,1])+j*omega*C_s)/(2*R_x)), np.sqrt((np.real(Y_far[i,0,1])+j*omega*C_s)/(2*R_x)))

##Run this only if P33 is needed  
Ct = W*N*epsilon_infty*1.207         
Ba_DART_RAM = np.imag(hilbert(Ga_DART_RAM))
Ba_DART_COM = np.imag(hilbert(Ga_DART_COM))
P_DART_RAM[:,2,2] = Ga_DART_RAM + j*(Ba_DART_RAM+omega*Ct)
P_DART_COM[:,2,2] = Ga_DART_COM + j*(Ba_DART_COM+omega*Ct)

###        P13 plotting
#        axs[i_1,i_2].plot(frequencies*1e-9,np.abs(P_DART_COM[:,0,2]), 'r', label='P13 DART - COM')
#        axs[i_1,i_2].plot(frequencies*1e-9,np.abs(P_DART_RAM[:,0,2]), 'b', label='P13 DART - RAM')
#        axs[i_1,i_2].plot(frequencies*1e-9,np.abs(P_double[:,0,2]), 'g', label='P13 double')
#        axs[i_1,i_2].legend()
#        axs[i_1,i_2].set_title(r'$f_c=$%.2f [GHz], N=%d' %(1e-9*v/lambda_0s[i_1],N))
#        fig.suptitle(r'P13 plots, R=%.4f j' %np.imag(R))
   
"""      S-matrix calculation and plotting """  
S_DART_RAM_iDART_toward = np.zeros((num_frequencies,2,2), dtype=np.complex)
S_DART_COM_iDART_toward = np.zeros((num_frequencies,2,2), dtype=np.complex)
S_DART_RAM_iDART_away = np.zeros((num_frequencies,2,2), dtype=np.complex)
S_DART_COM_iDART_away = np.zeros((num_frequencies,2,2), dtype=np.complex)

S_total_away = np.zeros((num_frequencies,2,2), dtype=np.complex)
S_total_toward = np.zeros((num_frequencies,2,2), dtype=np.complex)
R1=R2=Y0 = omega_c*W*C_S/K_squared
for i in range(num_frequencies):
#    P_iDART_RAM[i,:,:] = DART_to_iDART(P_DART_RAM[i,:,:])
    P_iDART_COM[i,:,:] = DART_to_iDART(P_DART_COM[i,:,:])
#    S_DART_RAM_iDART_toward[i,:,:] = Y2S_matrix(P2Y_matrix(P_DART_RAM[i,:,:], P_iDART_RAM[i,:,:]), R1, R2)
    S_DART_COM_iDART_toward[i,:,:] = Y2S_matrix(P2Y_matrix(P_DART_COM[i,:,:], P_iDART_COM[i,:,:]), R1, R2)
#    S_DART_RAM_iDART_away[i,:,:] = Y2S_matrix(P2Y_matrix(P_iDART_RAM[i,:,:], P_DART_RAM[i,:,:]), R1, R2)
    S_DART_COM_iDART_away[i,:,:] = Y2S_matrix(P2Y_matrix(P_iDART_COM[i,:,:], P_iDART_COM[i,:,:]), R1, R2)

    S_total_away[i,:,:] = Y2S_matrix(P2Y_matrix(P_DART_COM[i,:,:], P_DART_COM[i,:,:])+Z2Y_matrix(Z_far[i,:,:]), np.sqrt(((1/2*R_x)+j*omega*C_s)/(2*R_x)), np.sqrt(((1/2*R_x)+j*omega*C_s)/(2*R_x)))
    S_total_toward[i,:,:] = Y2S_matrix(P2Y_matrix(P_iDART_COM[i,:,:], P_iDART_COM[i,:,:])+Z2Y_matrix(Z_far[i,:,:]), np.sqrt(((1/2*R_x)+j*omega*C_s)/(2*R_x)), np.sqrt(((1/2*R_x)+j*omega*C_s)/(2*R_x))) 
    
fig, axs = plt.subplots(2,2)

###Broad frequancy range    
#axs[0,1].plot(frequencies*1e-9,20*np.log10(abs2(S_DART_COM_iDART_away[:,1,0]+S_far[:,1,0])), 'r', label='Away - COM')
#axs[0,1].plot(frequencies*1e-9,20*np.log10(abs2(S_DART_COM_iDART_toward[:,1,0]+S_far[:,1,0])), 'g', label='Toward DART - COM')
#axs[0,0].plot(frequencies*1e-9,20*np.log10(abs2(S_DART_COM_iDART_away[:,0,0]+S_far[:,0,0])), 'r', label='Away - COM')
#axs[0,0].plot(frequencies*1e-9,20*np.log10(abs2(S_DART_COM_iDART_toward[:,0,0]+S_far[:,0,0])), 'g', label='Toward - COM')

axs[0,1].plot(frequencies*1e-9,20*np.log10(abs2(S_total_away[:,1,0])), 'r', label='Away - COM')
axs[0,1].plot(frequencies*1e-9,20*np.log10(abs2(S_total_toward[:,1,0])), 'g', label='Toward DART - COM')
axs[0,0].plot(frequencies*1e-9,20*np.log10(abs2(S_total_away[:,0,0])), 'r', label='Away - COM')
axs[0,0].plot(frequencies*1e-9,20*np.log10(abs2(S_total_toward[:,0,0])), 'g', label='Toward - COM')

#ax.plot(frequencies*1e-9,20*np.log10(np.abs(S_DART_RAM_iDART_away[:,1,0])), 'k', label='Away DART - RAM')
#ax.plot(frequencies*1e-9,20*np.log10(np.abs(S_DART_RAM_iDART_toward[:,0,0])), 'b', label='Toward DART - RAM')
#ax.plot(frequencies*1e-9,20*np.log10(np.abs(S_DART_RAM_iDART_away[:,1,0])), 'k', label='Away DART - RAM')
#ax.plot(frequencies*1e-9,20*np.log10(np.abs(S_DART_RAM_iDART_toward[:,0,0])), 'b', label='Toward DART - RAM')

axs[0,0].legend()
axs[0,1].set_ylabel(r'$|S_{21}|^2$ [dB]')
axs[0,0].set_ylabel(r'$|S_{11}|^2$ [dB]')
axs[0,1].set_xlabel(r'Frequency [GHz]')
axs[0,0].set_xlabel(r'Frequency [GHz]')

###Near center frequency
#axs[1,1].plot(frequencies*1e-9,20*np.log10(abs2(S_DART_COM_iDART_away[:,1,0]+S_far[:,1,0])), 'r', label='Away - COM')
axs[1,1].plot(frequencies*1e-9,20*np.log10(abs2(S_total_away[:,1,0])), 'r', label='Away - COM')
#ax.plot(frequencies*1e-9,20*np.log10(np.abs(S_DART_RAM_iDART_away[:,1,0])), 'k', label='Away DART - RAM')
axs[1,1].plot(frequencies*1e-9,20*np.log10(abs2(S_total_toward[:,1,0])), 'g', label='Toward DART - COM')
#ax.plot(frequencies*1e-9,20*np.log10(np.abs(S_DART_RAM_iDART_toward[:,0,0])), 'b', label='Toward DART - RAM')
axs[1,0].plot(frequencies*1e-9,20*np.log10(abs2(S_total_away[:,0,0])), 'r', label='Away - COM')
#ax.plot(frequencies*1e-9,20*np.log10(np.abs(S_DART_RAM_iDART_away[:,1,0])), 'k', label='Away DART - RAM')
axs[1,0].plot(frequencies*1e-9,20*np.log10(abs2(S_total_toward[:,0,0])), 'g', label='Toward - COM')
#ax.plot(frequencies*1e-9,20*np.log10(np.abs(S_DART_RAM_iDART_toward[:,0,0])), 'b', label='Toward DART - RAM')
axs[1,1].set_ylabel(r'$|S_{21}|^2$ [dB]')
axs[1,0].set_ylabel(r'$|S_{11}|^2$ [dB]')
axs[1,1].set_xlabel(r'Frequency [GHz]')
axs[1,0].set_xlabel(r'Frequency [GHz]')
width_shown = 0.3e9 #[Hz]
axs[1,0].set_xlim([(f_c-width_shown/2)*1e-9, (f_c+width_shown/2)*1e-9])
axs[1,1].set_xlim([(f_c-width_shown/2)*1e-9, (f_c+width_shown/2)*1e-9])

fig.suptitle(r'S plots, R=%.4f j, $f_c=$%.2f [GHz], N=%d' %(np.imag(R), 1e-9*v/lambda_0,N))

###Directivity calculation and plotting
Dcom = np.abs(P_DART_COM[:,0,2]/P_DART_COM[:,1,2])
#Dcom = np.sqrt(np.abs(S_DART_COM_iDART_toward[:,1,0]/S_DART_COM_iDART_away[:,1,0]))


#fig2, axDcom = plt.subplots()
#axDcom.plot(frequencies*1e-9,20*np.log10(Dcom), '--', color='orange')
#axDcom.set_ylabel('Directivity - COM [dB]')
#axDcom.set_xlabel(r'Frequency [GHz]')
#fig2.suptitle(r'R=%.4f j, $f_c=$%.2f [GHz], N=%d' %(np.imag(R), 1e-9*v/lambda_0,N))

#        axs.set_xlabel('Frequency [GHz]')        

#plot of behavior far from center frequency (capacitive stray effect)
fig3, axs_far = plt.subplots(2)
axs_far[0].plot(frequencies*1e-9,20*np.log10(abs2(S_far[:,1,0])), '-')
axs_far[1].plot(frequencies*1e-9,20*np.log10(abs2(S_far[:,0,0])), '-')
axs_far[1].set_ylabel(r'$|S_{11}|^2$ [dB]')
axs_far[0].set_ylabel(r'$|Sfar_{21}|^2$ [dB]')
axs_far[1].set_xlabel(r'Frequency [GHz]')
axs_far[0].set_xlabel(r'Frequency [GHz]')
fig3.suptitle(r'R=%.4f j, $f_c=$%.2f [GHz], N=%d' %(np.imag(R), 1e-9*v/lambda_0,N))


plt.show()
print("--- %s seconds ---" % (time.time() - start_time))





   
    