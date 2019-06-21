# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:11:06 2019

@author: manip.batm
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy.special import legendre, lpmn
#from scipy.fftpack import fft, ifft

def legendre_function(s,x):
    a = 1
    Ps_x = 1
    for m in range(1,100):
        a = a*(m-1-s)*(m+s)*(1-x)/(2*m**2)
        Ps_x += a
    return Ps_x 

frequency = []
S11 = []
S12 = []
S21 = []

with open('C:/Users/manip.batm/Desktop/IDTs/Regular IDTs/From Jun-Liang/131.s2p', "r") as f:
    lines = f.readlines()
 
    # Loop through all lines, ignoring header.
    # Add last element to list (i.e. the process name)
    for l in lines[9:]:
        frequency.append(float(l.split()[0]))
        S11.append(float(l.split()[1]))
        S12.append(float(l.split()[3]))
        S21.append(float(l.split()[5]))

frequency = np.asarray(frequency)
S11 = np.asarray(S11)
S12 = np.asarray(S12)
S21 = np.asarray(S21)
#plt.plot(frequency,S11, label='S11')
#plt.plot(frequency,S12, label='S12')
#plt.plot(frequency,S21, label='S21')
#
FS12 = np.fft.rfft(S12)
n = len(frequency)
#freq_step = (frequency[-1]-frequency[0])/(n-1)
freq_step = (frequency[1]-frequency[0])
time = np.fft.rfftfreq(n, d=freq_step)
#plt.plot(time,np.abs(FS12), label='FS12')
##plt.xlim([5,95])

#Filterring

#tmin = 3e-7 #for 141
#tmax = 4.1e-7 #for 141
tmin = 3.4e-7 #for 131
tmax = 4.05e-7 #for 131
index= []
#for i,t in enumerate(time):
#    if t < tmin or t > tmax:
#        index.append(i)
#FS12filtered = np.delete(FS12,index) 
#timefiltered = np.delete(time,index) 
for i,t in enumerate(time):
    if t < tmin or t > tmax:
        FS12[i] = 0 
S12filtered = np.fft.irfft(FS12)
n = len(S12filtered)
time_step = (time[1]-time[0])
frequencyfiltered = np.fft.fftfreq(n, d=time_step)
plt.plot(frequency,np.abs(S12filtered), label='S12 filtered')
#plt.plot(timefiltered,np.abs(FS12filtered))

S12filtered = np.fft.irfft(FS12)
n = len(S12filtered)
time_step = (time[1]-time[0])
frequencyfiltered = np.fft.fftfreq(n, d=time_step)
#plt.plot(frequencyfiltered,np.abs(S12filtered), label='S12 filtered')
#plt.plot(frequencyfiltered,20*np.log10(np.abs(S12filtered)), label='S12 filtered')
plt.xlim(0.1,4.7*1e9)



#x= np.linspace(0,1000e-9,1000)
#y = np.sin(2*10e9*np.pi*x)
#xf =  np.fft.fftfreq(, d=timestep)
#yf = fft(y)
#plt.plot(x,yf)


#
#### Natural constants
#epsilon_0 = 8.8541878e-12
#j = np.complex(0,1)
#
####Design constants in S.I units
#lambda_0 = 700e-9/0.0158 #wavelength
#v = lambda_0 * 71e6 #wave velocity in the material
#W = 10*lambda_0 #aperture
#N = 80 #number of cells
#p = lambda_0 /4 #pitch
#eta = 5/8 #mettalization ratio
#a = p*eta #mettalized electrode width
#L = N*lambda_0 #total length
#omega_c = v*2*np.pi/lambda_0 #center frequency
#
#
#epsilon_infty = 3.8*epsilon_0 #effective permittivity
#Gamma_s = 0.01/epsilon_infty #coupling constant
#
#num_frequencies = 150
#frequencies = np.linspace(69e6,73e6,num_frequencies) #sent frequencies
#
#
#
#P = np.zeros((3,3), dtype=np.complex) #P-matrix
#c = np.zeros(N+1,dtype=np.complex) #forward wave
#b = np.zeros(N+1,dtype=np.complex) #backwards wave
#I = np.zeros(N+1,dtype=np.complex) #current
#P11 = np.zeros(num_frequencies)
#P12 = np.zeros(num_frequencies)
#Dram = np.zeros(num_frequencies)
#Dcom = np.zeros(num_frequencies)
#Ga = np.zeros(num_frequencies)
#c12 = np.zeros(num_frequencies, dtype=np.complex)
#conservation_check = np.zeros(num_frequencies, dtype=np.complex)
#
#R = np.complex(0,-0.0115) #reflection constant
#T = np.sqrt(1-np.abs(R)**2) #transmission constant
#
#def rho_f(k,a,p,epsilon_infty):
#    m = int((k*p)/(2*np.pi))
#    s = ((k*p)/(2*np.pi))-m
#    Delta = np.pi*a/p
#    return epsilon_infty * 2 *np.sin(np.pi*s)* (legendre_function(m,np.cos(Delta))/legendre_function(-s,-np.cos(Delta)))
#
####RAM equations for DART
#for i,f in enumerate(frequencies):
#    omega = f*2*np.pi
#    k = omega/v #wavenumber
#    
#    c[0] = 0
#    b[0] = 1
#    for n in range(1,N+1):
#        c[n] = (1/T)*c[n-1]*np.exp(-j*k*lambda_0) + (R/T)*b[n-1]
#        b[n] = -(R/T)*c[n-1] + (1/T)*b[n-1]*np.exp(j*k*lambda_0)
#        I[n] = -j*omega*W*rho_f(k,a,p,epsilon_infty)*(c[n]*np.exp(-j*k*(p/2))+b[n]*np.exp(j*k*(p/2)))*np.sqrt((2*Gamma_s)/(omega*W))
#    I_total = np.sum(I)
#        
#    P[1,1] = c[N]/b[N]
#    P[0,1] = b[0]/b[N]
#    P[2,1] = I_total/b[N]
#    P[1,2] = -P[2,1]/2
#    
#    c[N] = 1
#    b[N] = 0
#    for n in range(1,N+1):
#        c[-(n+1)] = (1/T)*c[-n]*np.exp(j*k*lambda_0) - (R/T)*b[-n]
#        b[-(n+1)] = (R/T)*c[-n] + (1/T)*b[-n]*np.exp(-j*k*lambda_0)
#        I[-(n+1)] = -j*omega*W*rho_f(k,a,p,epsilon_infty)*(c[-(n+1)]*np.exp(-j*k*(p/2))+b[-(n+1)]*np.exp(j*k*(p/2)))*np.sqrt((2*Gamma_s)/(omega*W))
#    I_total = np.sum(I)
#        
#    P[0,0] = b[0]/c[0]
#    P[1,0] = c[N]/c[0]
#    P[2,0] = I_total/c[0]
#    P[0,2] = -P[2,0]/2
#
#    Dram[i] = np.abs(P[0,2]/P[1,2])
#    Ga[i] = np.abs(P[0,2])**2 + np.abs(P[1,2])**2
#    P11[i] = np.abs(P[0,0])
#    P12[i] = np.abs(P[0,1])
#    conservation_check[i] = (P[0,0]*np.conj(P[0,2])+P[0,1]*np.conj(P[1,2]+P[0,2]))
# 
##    #for COM analysis
##    c12 = np.conj( -(R/lambda_0)*np.exp(-2*j*k*(lambda_0*(3/4))) )
##    aT = rho_f(k,a,p,epsilon_infty)*np.sqrt(omega*W*Gamma_s/2)
##    alpha1 = np.conj(-j*(aT/lambda_0)*np.exp(-j*k*lambda_0*(3/8)))
##    theta = np.angle(alpha1)
##    phi = np.angle(c12)
##    Dcom[i] = (np.abs(((1/np.tanh(np.abs(c12)*L/2))+np.exp(j*(2*theta-phi)))/((1/np.tanh(np.abs(c12)*L/2))-np.exp(-j*(2*theta-phi)))))**N
#
#
#Ga_c = omega_c*(epsilon_infty**2)*(N**2)*W*Gamma_s*1.556    
#plt.plot(frequencies*1e-6,20*np.log10(Dram), label='Directivity RAM')
##plt.plot(frequencies*1e-6,20*np.log10(Dcom), label='Directivity COM')
##plt.plot(frequencies*1e-6,Ga, label=r'Conductance $G_a$')
##plt.plot(omega_c/(2*np.pi)*1e-6,Ga_c, 'ro' , label=r'Center Conductance $G_a(f_c)$')
##plt.ylabel(r'$G_a$[mho]')
##plt.plot(frequencies*1e-6,np.abs(conservation_check), label=r'Conservation check')
#plt.xlabel('Frequency [MHz]')

    
plt.legend()
plt.show()