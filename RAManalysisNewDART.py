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
    for m in range(1,100):
        a = a*(m-1-s)*(m+s)*(1-x)/(2*m**2)
        Ps_x += a
    return Ps_x 

def rho_f(k,a,p,epsilon_infty):
    m = int((k*p)/(2*np.pi))
    s = ((k*p)/(2*np.pi))-m
    Delta = np.pi*a/p
    return epsilon_infty * 2 *np.sin(np.pi*s)* (legendre_function(m,np.cos(Delta))/legendre_function(-s,-np.cos(Delta)))

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

### Natural constants
epsilon_0 = 8.8541878e-12
j = np.complex(0,1)

###Material (GaAs) onstants
epsilon_infty = 10.89*epsilon_0 #effective permittivity
Gamma_s = 0.01/epsilon_infty #coupling constant
v = 2770 #wave velocity in the material


#lambda_0s = np.asarray([0.5,0.75,1,2])*1e-6 #wavelength
lambda_0s = np.asarray([1])*1e-6 #wavelength
#    L_times = np.array([20,40,60])*1e-9 #total length in s
L_times = np.array([40])*1e-9 #total length in s

Half_widths = np.zeros((len(lambda_0s),len(L_times)))
Max_amplitudes = np.zeros((len(lambda_0s),len(L_times)))
    
for i_1, lambda_0 in enumerate(lambda_0s):
    fig, ax = plt.subplots()

    Ls = v*L_times
    for i_2,L in enumerate(Ls):
        ###Design constants in S.I units 
        f_c = v/lambda_0 #center frequency 
        W = 30*1e-6 #aperture
        N = int(L/lambda_0) #number of cells
        p = lambda_0 /4 #pitch
        eta = 5/8 #mettalization ratio
        a = p*eta #mettalized electrode width
        omega_c = 2*np.pi*f_c #center frequency
        
        
        
        
        num_frequencies = 50
        frequencies = np.linspace(f_c-(num_frequencies/2)*1e5,f_c+(num_frequencies/2)*1e5,num_frequencies) #sent frequencies
        
        
        
        P = np.zeros((3,3), dtype=np.complex) #P-matrix
        c = np.zeros(N+1,dtype=np.complex) #forward wave
        b = np.zeros(N+1,dtype=np.complex) #backwards wave
        I = np.zeros(N+1,dtype=np.complex) #current
        P11 = np.zeros(num_frequencies, dtype=np.complex)
        P12 = np.zeros(num_frequencies, dtype=np.complex)
        P13 = np.zeros(num_frequencies, dtype=np.complex)
        P21 = np.zeros(num_frequencies, dtype=np.complex)
        P22 = np.zeros(num_frequencies, dtype=np.complex)
        P23 = np.zeros(num_frequencies, dtype=np.complex)
        P31 = np.zeros(num_frequencies, dtype=np.complex)
        P32 = np.zeros(num_frequencies, dtype=np.complex)
        P33 = np.zeros(num_frequencies, dtype=np.complex)
        Dram = np.zeros(num_frequencies)
        Dcom = np.zeros(num_frequencies)
        Ga = np.zeros(num_frequencies)
        c12 = np.zeros(num_frequencies, dtype=np.complex)
        conservation_check = np.zeros(num_frequencies, dtype=np.complex)
        
        R = np.complex(0,-0.0115) #reflection constant
        T = np.sqrt(1-np.abs(R)**2) #transmission constant
        
        ###RAM equations for DART
        for i,f in enumerate(frequencies):
            omega = f*2*np.pi
            k = omega/v #wavenumber
            
            c[0] = 0
            b[0] = 1
            for n in range(1,N+1):
                c[n] = (1/T)*c[n-1]*np.exp(-j*k*lambda_0) + (R/T)*b[n-1]
                b[n] = -(R/T)*c[n-1] + (1/T)*b[n-1]*np.exp(j*k*lambda_0)
                I[n] = -j*omega*W*rho_f(k,a,p,epsilon_infty)*(c[n]*np.exp(-j*k*(p/2))+b[n]*np.exp(j*k*(p/2)))*np.sqrt((2*Gamma_s)/(omega*W))
            I_total = np.sum(I)
                
            P22[i] = P[1,1] = c[N]/b[N]
            P12[i] = P[0,1] = b[0]/b[N]
            P32[i] = P[2,1] = I_total/b[N]
            P23[i] = P[1,2] = -P[2,1]/2
            
            c[N] = 1
            b[N] = 0
            for n in range(1,N+1):
                c[-(n+1)] = (1/T)*c[-n]*np.exp(j*k*lambda_0) - (R/T)*b[-n]
                b[-(n+1)] = (R/T)*c[-n] + (1/T)*b[-n]*np.exp(-j*k*lambda_0)
                I[-(n+1)] = -j*omega*W*rho_f(k,a,p,epsilon_infty)*(c[-(n+1)]*np.exp(-j*k*(p/2))+b[-(n+1)]*np.exp(j*k*(p/2)))*np.sqrt((2*Gamma_s)/(omega*W))
            I_total = np.sum(I)
                
            P11[i] = P[0,0] = b[0]/c[0]
            P21[i] = P[1,0] = c[N]/c[0]
            P31[i] = P[2,0] = I_total/c[0]
            P13[i] = P[0,2] = -P[2,0]/2
    
            Dram[i] = np.abs(P[0,2]/P[1,2])
            Ga[i] = np.abs(P[0,2])**2 + np.abs(P[1,2])**2
#           conservation_check[i] = (P[0,0]*np.conj(P[0,2])+P[0,1]*np.conj(P[1,2]+P[0,2]))
         
        #    #for COM analysis
        #    c12 = np.conj( -(R/lambda_0)*np.exp(-2*j*k*(lambda_0*(3/4))) )
        #    aT = rho_f(k,a,p,epsilon_infty)*np.sqrt(omega*W*Gamma_s/2)
        #    alpha1 = np.conj(-j*(aT/lambda_0)*np.exp(-j*k*lambda_0*(3/8)))
        #    theta = np.angle(alpha1)
        #    phi = np.angle(c12)
        #    Dcom[i] = (np.abs(((1/np.tanh(np.abs(c12)*L/2))+np.exp(j*(2*theta-phi)))/((1/np.tanh(np.abs(c12)*L/2))-np.exp(-j*(2*theta-phi)))))**N

#Run this only if P33 is needed           
        Ba = np.imag(hilbert(Ga))
        Ct = W*N*epsilon_infty*1.207
        P33 = Ga + j*(Ba+omega*Ct) 

        Ga_c = omega_c*(epsilon_infty**2)*(N**2)*W*Gamma_s*1.556    
        ax.plot(frequencies*1e-9,20*np.log10(Dram), label=r'N=%d, $f_c$=%.2f[Ghz]' %(N,f_c*1e-9 ))
        Max_amplitudes[i_1,i_2] = max(20*np.log10(Dram))
        Half_widths[i_1,i_2] = np.array(FWHM(np.array(frequencies*1e-9),np.array(20*np.log10(Dram))))
        #plt.plot(frequencies*1e-9,20*np.log10(Dcom), label='Directivity COM')
        #plt.plot(frequencies*1e-9,Ga, label=r'Conductance $G_a$')
        #plt.plot(omega_c/(2*np.pi)*1e-6,Ga_c, 'ro' , label=r'Center Conductance $G_a(f_c)$')
        #plt.ylabel(r'$G_a$[mho]')
        #plt.plot(frequencies*1e-6,np.abs(conservation_check), label=r'Conservation check')
    ax.set_xlabel('Frequency [GHz]')
    
    ax.set_title('Directivity RAM')
        
    ax.legend()
    plt.show()
    print("--- %s seconds ---" % (time.time() - start_time))

####Plot amplitudes and bandwidths as a function of DART length and center frequency
#for i in range(lambda_0s.size):    
#    fig, axs = plt.subplots(2,1)
#    fig.suptitle(r'Constant $f_c=$%.2f [GHz]' %(1e-9*v/lambda_0s[i]))       
#    axs[0].plot(L_times*1e9, Max_amplitudes[i,:], 'bo')
#    axs[0].set_ylabel('Max Directivity [dB]')
#    axs[1].plot(L_times*1e9, Half_widths[i,:], 'bo')
#    axs[1].set_ylabel('Band Width [GHz]')
#    axs[1].set_xlabel('DART length [ns]')
#           
#    plt.show()
#
#for j in range(L_times.size):    
#    fig, axs = plt.subplots(2,1)
#    fig.suptitle(r'Constant $L=$%.f [ns]' %(L_times[j]*1e9))       
#    axs[0].plot(1e-9*v/lambda_0s, Max_amplitudes[:,j], 'bo')
#    axs[0].set_ylabel('Max Directivity [dB]')
#    axs[1].plot(1e-9*v/lambda_0s, Half_widths[:,j], 'bo')
#    axs[1].set_ylabel('Band Width [GHz]')
#    axs[1].set_xlabel('$f_c$ [GHz]')
#           
#    plt.show()



   
    