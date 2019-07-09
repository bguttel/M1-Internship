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

def abs2(x):
    return x.real**2 + x.imag**2

def legendre_function(s,x):
    a = 1
    Ps_x = 1
    for m in range(1,100):
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

#def Q_m_integrand(s,m,a,p):
#    """Gives the chrage on the m_th electrode"""
#    Delta = np.pi*a/p
#    integrand = np.sin(np.pi*s)*np.cos(2*np.pi*m*s)*(legendre_function(-s,np.cos(Delta))/legendre_function(-s,-np.cos(Delta)))
#    return integrand
#
#def integrate(m, a, p, init, end, dx=0.1):
#    i = init
#    integral = 0
#    while i <= end:
#        integral += Q_m_integrand(i,m,a,p)*dx
#        i += dx
#    return integral
#
#def C_t(W,N,epsilon_infty,a,p):
#    """"Gives the total capacitance of a DART"""
#    C = 0
#    P = np.zeros(N)
#    for i in range(N): #DART electrode configuration
#        if i % 3 == 0:
#            P[i] =1
#    for i in range(N):
#        for j in range(N):
#            C += P[i]*P[j]*integrate(np.abs(i-j),a,p,init=0,end=1, dx=0.1)
#    return 2*epsilon_infty*W*C

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

def P_matrix_double(k,lambda_0,a,p,epsilon_infty,omega,W,Gamma_s,N):
    """P-matrix of double-fingered IDT (except P33). Based on Morgan's book, p.132"""
    L = N*lambda_0
    S_e = 4 #4 electrodes per IDT cell
    P = np.zeros((3,3), dtype=np.complex) #P-matrix
    P[1,0] = P[0,1] = np.exp(-j*k*L)
    P[0,2] = j*rho_e(N,S_e,k,a,p,epsilon_infty)*np.sqrt(omega*W*(Gamma_s/2))*np.exp(-j*k*L/2)
    P[1,2] = j*rho_e(N,S_e,-k,a,p,epsilon_infty)*np.sqrt(omega*W*(Gamma_s/2))*np.exp(-j*k*L/2)
    P[2,0] = -2*P[0,2]
    P[2,1] = -2*P[1,2]
    return P

def P_matrix_double2(k,lambda_0,a,p,epsilon_infty,omega,W,Gamma_s,N,C_S,omega_c,K_squared):
    """P-matrix of double-fingered IDT. Based on Aref's paper, https://arxiv.org/pdf/1506.01631.pdf"""
    L = N*lambda_0
    S_e = 4 #4 electrodes per IDT cell
    P = np.zeros((3,3), dtype=np.complex) #P-matrix
    P[1,0] = P[0,1] = np.exp(-j*k*L)
    P[0,2] = j*rho_e(N,S_e,k,a,p,epsilon_infty)*np.sqrt(omega*W*(Gamma_s/2))*np.exp(-j*k*L/2)
    P[1,2] = j*rho_e(N,S_e,-k,a,p,epsilon_infty)*np.sqrt(omega*W*(Gamma_s/2))*np.exp(-j*k*L/2)
    P[2,0] = -2*P[0,2]
    P[2,1] = -2*P[1,2]
    
    Ct = np.sqrt(2)*N*C_S*W
    X = N*np.pi*(omega-omega_c)/omega_c
    Ga0 = 1.3*K_squared*(N**2)*omega_c*W*C_S
    Ga = Ga0*np.square(np.sin(X)/X)
    Ba = Ga0*(np.sin(2*X)-2*X)/2*np.square(X)
    P[2,2] = Ga +j*(Ba+omega*Ct)
    
    return P

def P2Y_matrix(Pa,Pb):
    """Y-matrix of a delay-line using P-matrices"""
    Y = np.zeros((2,2), dtype=np.complex) #P-matrix
    Y[0,0] = Pa[2,2] - 2*Pb[0,0]*np.square(Pa[0,2])/(1-Pa[0,0]*Pb[0,0])
    Y[1,1] = Pb[2,2] - 2*Pa[0,0]*np.square(Pb[0,2])/(1-Pa[0,0]*Pb[0,0])
    Y[0,1] = Y[1,0] = -2*Pa[0,2]*Pb[0,2]/(1-Pa[0,0]*Pb[0,0])
    
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
    
    
    
    
### Natural constants
epsilon_0 = 8.8541878e-12 #[F/m]
j = np.complex(0,1)

###Material (GaAs) constants
epsilon_infty = 10.89*epsilon_0 #effective permittivity
K_squared = 7e-4 #Piezo-electric coupling coefficient
Gamma_s = (K_squared/2)/epsilon_infty #coupling constant
v_f = 2770 #free wave velocity in the material
v = v_f #wave velocity in the material


C_S = 1.2e-10 #[F/m] Effective capacitance 


lambda_0s = np.asarray([0.5,0.75,1,1.5,2])*1e-6 #wavelength
#lambda_0s = np.asarray([0.5])*1e-6 #wavelength
Ns = np.array([90,100,110,120,130,140,150]) #total number of periods
#Ns = np.array([130]) #total number of periods
#Rs = np.array([-0.05,-0.02,0.02,0.05])*j
Rs = np.array([0.02])*j

num_frequencies = 200

            

Half_widths = np.zeros((lambda_0s.size,Ns.size))
Max_amplitudes = np.zeros((lambda_0s.size,Ns.size))


for R in Rs: #reflection constant
    T = np.sqrt(1-abs2(R)) #transmission constant
#    fig, axs = plt.subplots(lambda_0s.size,Ns.size)    
    for i_1, lambda_0 in enumerate(lambda_0s):
        
    
        for i_2,N in enumerate(Ns):
            ###Design constants in S.I units 
            f_c = v/lambda_0 #center frequency 
            W = 30*1e-6 #aperture
            L = N*lambda_0 #total langth
            p = lambda_0 /4 #pitch
            eta = 0.5 #mettalization ratio
            a = p*eta #mettalized electrode width
            omega_c = 2*np.pi*f_c #center frequency
            k_c = omega_c / v
            
            
            frequencies_width = 1 * 1e9 #Desired total width of frequencies simulated, <2*f_c
            if frequencies_width > 2*f_c:
                print(r'frequencies width  = %.2f > 2* %.2f = $f_c$' %(frequencies_width*1e-9, f_c*1e-9))
                frequencies_width = 2*f_c - 1
            frequencies = np.linspace(f_c-(frequencies_width/2),f_c+(frequencies_width/2),num_frequencies) #sent frequencies
            
            
            
            """P-matrix calculations using RAM and COM analyses"""
            P_DART_RAM = np.zeros((num_frequencies,3,3), dtype=np.complex) #RAM DART P-matrix
            P_DART_COM = np.zeros((num_frequencies,3,3), dtype=np.complex) #COM DART P-matrix
            P_double = np.zeros((num_frequencies,3,3), dtype=np.complex) #double IDT P-matrix
            c = np.zeros(N+1, dtype=np.complex) #forward wave
            b = np.zeros(N+1, dtype=np.complex) #backwards wave
            I = np.zeros(N+1, dtype=np.complex) #current
            
            P13_double = np.zeros(num_frequencies, dtype=np.complex)
            P33_double = np.zeros(num_frequencies, dtype=np.complex)
            Dram = np.zeros(num_frequencies)
            Dcom = np.zeros(num_frequencies)
            Ga_DART_RAM = np.zeros(num_frequencies)
            Ga_DART_COM = np.zeros(num_frequencies)
            Ga_double = np.zeros(num_frequencies)
    #        c12 = np.zeros(num_frequencies, dtype=np.complex)
            Dram = np.zeros(num_frequencies)
            Dcom = np.zeros(num_frequencies)
            Ga = np.zeros(num_frequencies)
            conservation_check = np.zeros(num_frequencies, dtype=np.complex)
            
    ###     RAM analysis
            for i,f in enumerate(frequencies):
                omega = f*2*np.pi
 ###             Velocity adjustment (elctrical loading)               
                v_m = (1-K_squared/2)*v_f #[m/s] velocity on metal
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
        
                Dram[i]= np.abs(P_DART_RAM[i,1,2]/P_DART_RAM[i,0,2])
                Ga_DART_RAM[i] = np.abs(P_DART_RAM[i,0,2])**2 + np.abs(P_DART_RAM[i,1,2])**2
    #           conservation_check[i] = (P[0,0]*np.conj(P[0,2])+P[0,1]*np.conj(P[1,2]+P[0,2]))
             
    ###         COM analysis
                c12 = np.conj( -(R/lambda_0)*np.exp(-2*j*k*(lambda_0*(3/4))) )
                aT = rho_f(k,a,p,epsilon_infty)*np.sqrt(omega*W*Gamma_s/2)
                alpha1 = np.conj(-j*(aT/lambda_0)*np.exp(-j*k*lambda_0*(3/8)))
                theta = np.angle(alpha1)
                phi = np.angle(c12)
                P_DART_COM[i,0,0] = -(np.conj(c12)/np.abs(c12))*np.tanh(np.abs(c12)*L)
                P_DART_COM[i,0,1] = np.exp(-j*k*L)/np.cosh(np.abs(c12)*L)
                delta = k - k_c
                s = np.sqrt(np.complex(delta**2 - np.abs(c12)**2))
                K1 = (np.conj(alpha1)*c12-j*delta*alpha1)/(s**2)
                K2 = (np.conj(c12)*alpha1+j*delta*np.conj(alpha1))/(s**2)
                D = s*np.cos(s*L) + j*delta*np.sin(s*L)
                P_DART_COM[i ,2,0] = (2*np.conj(alpha1)*np.sin(s*L) - 2*s*K2*(np.cos(s*L)-1))/D
                P_DART_COM[i,2,1] = np.exp(-j*k_c*L)*((-2*alpha1*np.sin(s*L) - 2*s*K1*(np.cos(s*L)-1))/D)
                P_DART_COM[i,0,2] = P_DART_COM[i,2,0]/(-2)
                P_DART_COM[i,1,2] = P_DART_COM[i,2,1]/(-2)
                Ga_DART_COM[i] = np.abs(P_DART_COM[i,0,2])**2 + np.abs(P_DART_COM[i,1,2])**2
#                Dcom[i] = (np.abs(((1/np.tanh(np.abs(c12)*L/2))+np.exp(j*(2*theta-phi)))/((1/np.tanh(np.abs(c12)*L/2))-np.exp(-j*(2*theta-phi)))))**N
                Dcom[i]= np.abs(P_DART_COM[i,1,2]/P_DART_COM[i,0,2])
    
    ###            double-fingered IDT P-matrix, two methods
    #            P_double[i,:,:] = P_matrix_double(k,lambda_0,a,p,epsilon_infty,omega,W,Gamma_s,N)
    #            Ga_double[i] = omega*(epsilon_infty**2)*(N**2)*W*Gamma_s*3.111
    #        Ct_double = np.ones(num_frequencies)*W*N*epsilon_infty*1.414
    #        P33_double = Ga_double+j*(np.imag(hilbert(Ga_double))+omega*Ct_double)
    #
    #        
    #        P_double[:,2,2] = Ga_double+j*(np.imag(hilbert(Ga_double))+omega*Ct_double)
            
#                P_double[i,:,:] = P_matrix_double2(k,lambda_0,a,p,epsilon_infty,omega,W,Gamma_s,N,C_S,omega_c,K_squared)
            
            
    ##Run this only if P33 is needed  
            Ct = W*N*epsilon_infty*1.207         
            Ba_DART_RAM = np.imag(hilbert(Ga_DART_RAM))
            Ba_DART_COM = np.imag(hilbert(Ga_DART_COM))
            P_DART_RAM[:,2,2] = Ga_DART_RAM + j*(Ba_DART_RAM+omega*Ct)
            P_DART_COM[:,2,2] = Ga_DART_COM + j*(Ba_DART_COM+omega*Ct)
    
    #        Ga_c = omega_c*(epsilon_infty**2)*(N**2)*W*Gamma_s*1.556    
    #        ax.plot(frequencies*1e-9,20*np.log10(Dram), label=r'N=%d, $f_c$=%.2f[Ghz]' %(N,f_c*1e-9 ))
            Max_amplitudes[i_1,i_2] = max(20*np.log10(np.abs(Dcom)))
            Half_widths[i_1,i_2] = np.array(FWHM(np.array(frequencies*1e-9),np.array(20*np.log10(np.abs(Dcom)))))
            #plt.plot(frequencies*1e-9,20*np.log10(Dcom), label='Directivity COM')
            #plt.plot(frequencies*1e-9,Ga, label=r'Conductance $G_a$')
            #plt.plot(omega_c/(2*np.pi)*1e-6,Ga_c, 'ro' , label=r'Center Conductance $G_a(f_c)$')
            #plt.ylabel(r'$G_a$[mho]')
            #plt.plot(frequencies*1e-6,np.abs(conservation_check), label=r'Conservation check')
            
            """P13 plotting"""
#            axs[i_1,i_2].plot(frequencies*1e-9,np.abs(P_DART_COM[:,0,2]), 'r', label='P13 DART - COM')
#            axs[i_1,i_2].plot(frequencies*1e-9,np.abs(P_DART_RAM[:,0,2]), 'b', label='P13 DART - RAM')
#            axs[i_1,i_2].plot(frequencies*1e-9,np.abs(P_double[:,0,2]), 'g', label='P13 double')
#            axs[i_1,i_2].legend()
#            axs[i_1,i_2].set_title(r'$f_c=$%.2f [GHz], N=%d' %(1e-9*v/lambda_0s[i_1],N))
#            fig.suptitle(r'P13 plots, R=%.4f j' %np.imag(R))
  
            """"Directivity plotting"""
            fig, axs = plt.subplots()
            axs.plot(frequencies*1e-9,20*np.log10(np.abs(P_DART_COM[:,1,2]/P_DART_COM[:,0,2])), 'r', label='Directivity - COM')
#            axs.plot(frequencies*1e-9,20*np.log10(np.abs(P_DART_RAM[:,1,2]/P_DART_RAM[:,0,2])), '--b', label='Directivity - RAM')
#            axs.legend()
            axs.set_title(r'$f_c=$%.2f [GHz], N=%d, R=%.3fj ' %(1e-9*v_f/lambda_0s[i_1],N,np.imag(R)))
            axs.set_xlabel(r'Frequency [GHz]')
            axs.set_ylabel(r'Directivity [dB]')
            fig.suptitle(r'Directivity - COM')
            fig.savefig('Simulations/Directivities/Directivity_f%3.f_N%d' %(1e-7*v_f/lambda_0s[i_1],N))
            
            
            
            
            
            """S21 calculation and plotting"""
#            S21_DART_RAM_double = np.zeros(num_frequencies, dtype=np.complex)
#            S21_DART_COM_double = np.zeros(num_frequencies, dtype=np.complex)
#            S21_double_double = np.zeros(num_frequencies, dtype=np.complex)
#            R1=R2=Y0 = omega_c*W*C_S/K_squared
#            for i in range(num_frequencies):
#    #            S21_DART_RAM_double[i] = S21(P_DART_RAM[i,:,:],P_double[i,:,:], Y0)
#    #            S21_DART_COM_double[i] = S21(P_DART_COM[i,:,:],P_double[i,:,:], Y0)
#    #            S21_double_double[i] = S21(P_double[i,:,:],P_double[i,:,:], Y0)
#                
#                S21_DART_RAM_double[i] = Y2S_matrix(P2Y_matrix(P_DART_RAM[i,:,:],P_double[i,:,:]), R1, R2)[1,0]
#                S21_DART_COM_double[i] = Y2S_matrix(P2Y_matrix(P_DART_COM[i,:,:],P_double[i,:,:]), R1, R2)[1,0]
#                S21_double_double[i] = Y2S_matrix(P2Y_matrix(P_double[i,:,:],P_double[i,:,:]), R1, R2)[1,0]
#                
#            axs[i_1,i_2].plot(frequencies*1e-9,20*np.log10(abs2(S21_DART_RAM_double)), 'r', label='S21 DART - RAM')
#            axs[i_1,i_2].plot(frequencies*1e-9,20*np.log10(abs2(S21_DART_COM_double)), 'b', label='S21 DART - COM')
#            axs[i_1,i_2].plot(frequencies*1e-9,20*np.log10(abs2(S21_double_double)), 'g', label='S21 double')
#            axs[i_1,i_2].legend()
#            axs[i_1,i_2].set_title(r'$f_c=$%.2f [GHz], N=%d' %(1e-9*v_f/lambda_0s[i_1],N))
#            fig.suptitle(r'S21 plots, R=%.4f j' %np.imag(R))
    
    
        
#        plt.show()
        print("--- %s seconds ---" % (time.time() - start_time))
    
    """Plot amplitudes and bandwidths as a function of DART length and center frequency"""
    for i_1 in range(lambda_0s.size):    
        fig, axs = plt.subplots(2,1)
        fig.suptitle(r'Constant $f_c=$%.2f [GHz], R=%.3fj ' %(1e-9*v_f/lambda_0s[i_1], np.imag(R)))       
        axs[0].plot(Ns, Max_amplitudes[i_1,:], 'bo')
        axs[0].set_ylabel('Max Directivity [dB]')
        axs[1].plot(Ns, Half_widths[i_1,:], 'go')
        axs[1].set_ylabel(' Band Width [GHz]')
        axs[1].set_xlabel('Number of DART periods')
        fig.savefig('Simulations/MaxDirectivities_BandWidths/MaxDirectivities_BandWidths_f%3.f' %(1e-7*v_f/lambda_0s[i_1]))
#        fig.savefig(r'f_c%.2f[GHz]R%.3fj' %(1e-9*v/lambda_0s[i_1], np.imag(R)))
               
        plt.show()
    
    for i_2 in range(Ns.size):    
        fig, axs = plt.subplots(2,1)
        fig.suptitle(r'Constant N=%.f, R=%.3fj' %(Ns[i_2], np.imag(R)))       
        axs[0].plot(1e-9*v/lambda_0s, Max_amplitudes[:,i_2], 'bo')
        axs[0].set_ylabel('Max Directivity [dB]')
        axs[1].plot(1e-9*v/lambda_0s, Half_widths[:,i_2], 'go')
        axs[1].set_ylabel('Band Width [GHz]')
        axs[1].set_xlabel('$f_c$ [GHz]')
        fig.savefig('Simulations/MaxDirectivities_BandWidths/MaxDirectivities_BandWidths_N%d' %(Ns[i_2]))
               
        plt.show()



   
    