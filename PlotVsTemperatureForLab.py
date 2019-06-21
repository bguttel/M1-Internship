import numpy as np
import h5py
import matplotlib.pyplot as plt
from os.path import isfile, join
import os
import datetime as dt

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
            r -=1
    b = np.append(b,a[-1])
    return b

# Changables
data_path = 'C:/Users/manip.batm/Desktop/DClines/cd3Data'
filename_dummy = 'cd3_25-05-2019_15_20_43_1.h5'
init_date = dt.date(2019,6,4)
end_date = dt.date(2019,6,4)
init_time = dt.time(17,11,00)
end_time = dt.time(21,59,00)



init_combine = dt.datetime.combine(init_date,init_time)
end_combine = dt.datetime.combine(end_date,end_time)

T1 = []
T2 = []
T3 = []
T = []

V0 = []
V1 = []
V2 = []
V3 = []

R0 = []
R1 = []
R2 = []
R3 = []

adc0 = []
adc1 = []
adc2 = []
adc3 = []

adcs = [adc0,adc1,adc2,adc3]
Rs = [R0,R1,R2,R3]
Vs = [V0,V1,V2,V3]

adcs_name = ['ADC0','ADC1','ADC2','ADC3']
gates = ['0-1','0-2','0-3','0-4']

with open('C:/Users/manip.batm/Desktop/DClines/TemperatureData/thermo_2019-06-04.txt', "r") as f:
    lines = f.readlines()
 
    # Loop through all lines, ignoring header.
    # Add last element to list (i.e. the process name)
    for l in lines[1:]:
        TempDate = dt.datetime.strptime(l.split()[0],"%d/%m/%Y").date()
        TempTime = dt.datetime.strptime(l.split()[1],"%H:%M:%S").time()
        TempCombine = dt.datetime.combine(TempDate,TempTime)
        if   TempCombine >= init_combine and  TempCombine <= end_combine:
            T1.append(float(l.split()[6]))
        
with open('C:/Users/manip.batm/Desktop/DClines/TemperatureData/thermo_2019-06-05.txt', "r") as f:
    lines = f.readlines()
 
    # Loop through all lines, ignoring header.
    # Add last element to list (i.e. the process name)
    for l in lines[1:]:
        TempDate = dt.datetime.strptime(l.split()[0],"%d/%m/%Y").date()
        TempTime = dt.datetime.strptime(l.split()[1],"%H:%M:%S").time()
        TempCombine = dt.datetime.combine(TempDate,TempTime)
        if   TempCombine >= init_combine and  TempCombine <= end_combine:
            T1.append(float(l.split()[6]))

#with open('C:/Users/manip.batm/Desktop/27.5/DClines_cd1/TRMC_data/thermo_2019-05-28.txt', "r") as f:
#    lines = f.readlines()
# 
#    # Loop through all lines, ignoring header.
#    # Add last element to list (i.e. the process name)
#    for l in lines[1:]:
#        if dt.datetime.strptime(l.split()[1],"%H:%M:%S").time() <= end_time:
#            T3.append(float(l.split()[6]))
        
for root, dirs, files in os.walk(data_path):
    for file in files:
        if len(file) == len(filename_dummy):
            full_path = join(root,file)
            file = h5py.File(full_path,'r')
            for i,adc in enumerate(adcs):
                V = np.array(file['DRIVERS']['labview'][gates[i]]['values']).flatten()
                Vs[i].append(np.copy(V))
                d0 = file['data'][adcs_name[i]]
                d0 = np.einsum('ij->ji',d0).flatten()
                Rs[i].append(np.copy(abs(V/d0)))
                adc.append(np.copy(d0))
            file.close()


adc0 = np.array(adcs[0]).flatten()
adc1 = np.array(adcs[1]).flatten()
adc2 = np.array(adcs[2]).flatten()
adc3 = np.array(adcs[3]).flatten()

R0 = np.array(Rs[0]).flatten()
R1 = np.array(Rs[1]).flatten()
R2 = np.array(Rs[2]).flatten()
R3 = np.array(Rs[3]).flatten()

V0 = np.array(Vs[0]).flatten()
V1 = np.array(Vs[1]).flatten()
V2 = np.array(Vs[2]).flatten()
V3 = np.array(Vs[3]).flatten()
 

T1 = np.asarray(T1)
T2 = np.asarray(T2)
T3 = np.asarray(T3)
T = np.concatenate((T1,T2,T3), axis=None)
T = Stretch(T,len(adc0))


fig,ax = plt.subplots(1)

ax.plot(T,adc0,label='adc0')
ax.plot(T,adc1,label='adc1')
ax.plot(T,adc2,label='adc2')
ax.plot(T,adc3,label='adc3')
ax.set_ylabel(r'$ V_{out} \ ( V)$')
#
#ax.plot(T,R0/100,label='R0')
#ax.plot(T,R1/1000,label='R1')
#ax.plot(T,R2/1000,label='R2')
#ax.plot(T,R3/1000,label='R3')
#ax.set_ylabel(r'R $(\Omega)$')
#
ax.set_xlabel('tenperature (K)')
#ax.set_xlabel('time (min)')

ax.set_title('CoolDown_%s-%s' %(init_combine.strftime('%d/%m,%H:%M'), end_combine.strftime('%d/%m,%H:%M')))
# ax.set_ylim((0,2000))
ax.legend(loc='right')
fig.tight_layout()
plt.show()