import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure
import matplotlib.ticker as mtick
import os

pwd = os.getcwd()
pwd = pwd.split('/')

figure(figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')

X1, Y1 = [], []
X2, Y2 = [], []
data = []

data = np.loadtxt('ErrorMap.dat')

rho_err = data[:,  0]  
Ux_err  = data[:,  1]  
Uy_err  = data[:,  2]  
Uz_err  = data[:,  3]  
P_err   = data[:,  4]  
T_err   = data[:,  5]  
increment = data[:,  6]  

rho     = data[0,7]
Ux      = data[0,8]
Uy      = data[0,9]
Uz      = data[0,10]
LorentzFactor = data[0,12]
T       = data[:, 13]  

plt.plot(T, rho_err,  'bo', label='rho'  )
plt.plot(T, Ux_err,   'g^', label='Ux'   )
plt.plot(T, Uy_err,   'b<', label='Uy'   )
plt.plot(T, Uz_err,   'rx', label='Uz'   )
plt.plot(T, P_err,    'rD', label='P'    )
plt.plot(T, T_err,    'k+', label='T'  )
plt.xlabel('temperature (GeV)', size='25')
plt.ylabel('relative error'   , size='25')

x_value = T[np.where(np.isnan(Ux_err))[0]]

#ifdef  FLOAT8
plt.hlines(+np.finfo(float).eps,min(T),max(T),colors='r', linestyle='--')
plt.hlines(-np.finfo(float).eps,min(T),max(T),colors='r', linestyle='--')
plt.yscale('symlog', linthreshy=2e-16)
#else
#plt.hlines(+1.192092896e-07,min(T),max(T),colors='r', linestyle='--')
#plt.hlines(-1.192092896e-07,min(T),max(T),colors='r', linestyle='--')
#plt.yscale('symlog', linthreshy=1e-8)
#endif

#plt.vlines(x_value,-1,1,colors='b')

plt.tick_params(axis='both', labelsize='15')
plt.legend(fontsize='20', loc='upper right')
plt.ylim(-1e15,1e15)
#plt.xlim(1e-6,1e3)
plt.xscale('log')
#plt.yticks([-1e-16, -1e-15,0, 1e-16, 1e-15])
plt.title("%s (rho=%3.2e, Ux=%3.2e, Uy=%3.2e, Uz=%3.2e, LorentzFactor=%3.2e)" % (pwd[-1], rho, Ux, Uy, Uz, LorentzFactor), size='25')
plt.grid(True)
plt.show()
#plt.savefig('fun_1.png')
