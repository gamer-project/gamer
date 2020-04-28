import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import os

pwd = os.getcwd()
pwd = pwd.split('/')

f, ax = plt.subplots( 6, 2, sharex=False, sharey=False )
f.subplots_adjust( hspace=0.05, wspace=0.15 )
f.set_size_inches( 12.0, 10.0 )



X1, Y1 = [], []
X2, Y2 = [], []
data = []

data1 = np.loadtxt('ErrorMap1.dat')
data2 = np.loadtxt('ErrorMap2.dat')

rho1_err       = data1[:,0]  
Ux1_err        = data1[:,1]  
Uy1_err        = data1[:,2]  
Uz1_err        = data1[:,3]  
P1_err         = data1[:,4]  
T1_err         = data1[:,5]  
increment      = data1[:,6]  

rho1           = data1[0,7]
Ux1            = data1[0,8]
Uy1            = data1[0,9]
Uz1            = data1[0,10]
LorentzFactor1 = data1[0,12]
T1             = data1[:,13]  

rho2_err       = data2[:,0]  
Ux2_err        = data2[:,1]  
Uy2_err        = data2[:,2]  
Uz2_err        = data2[:,3]  
P2_err         = data2[:,4]  
T2_err         = data2[:,5]  
increment      = data2[:,6]  

rho2           = data2[0,7]
Ux2            = data2[0,8]
Uy2            = data2[0,9]
Uz2            = data2[0,10]
LorentzFactor2 = data2[0,12]
T2             = data2[:,13]  



ax[0][0].plot(T1, abs(T1_err  ),  '+', label=r'$\epsilon_{T}$'    , ms=4 , markerfacecolor='none', color='r' )
ax[1][0].plot(T1, abs(P1_err  ),  'D', label=r'$\epsilon_{P}$'    , ms=4 , markerfacecolor='none', color='y' )
ax[2][0].plot(T1, abs(rho1_err),  'o', label=r'$\epsilon_{\rho}$' , ms=4 , markerfacecolor='none', color='g' )
ax[3][0].plot(T1, abs(Ux1_err ),  '^', label=r'$\epsilon_{U_x}$'  , ms=4 , markerfacecolor='none', color='c' )
ax[4][0].plot(T1, abs(Uy1_err ),  '<', label=r'$\epsilon_{U_y}$'  , ms=4 , markerfacecolor='none', color='b' )
ax[5][0].plot(T1, abs(Uz1_err ),  'x', label=r'$\epsilon_{U_z}$'  , ms=4 , markerfacecolor='none', color='m' )


ax[0][1].plot(T2, abs(T2_err  )/1e-16,  '+', label=r'$\epsilon_{T}$'    , ms=4 , markerfacecolor='none', color='r' )
ax[1][1].plot(T2, abs(P2_err  )/1e-16,  'D', label=r'$\epsilon_{P}$'    , ms=4 , markerfacecolor='none', color='y' )
ax[2][1].plot(T2, abs(rho2_err)/1e-16,  'o', label=r'$\epsilon_{\rho}$' , ms=4 , markerfacecolor='none', color='g' )
ax[3][1].plot(T2, abs(Ux2_err )/1e-16,  '^', label=r'$\epsilon_{U_x}$'  , ms=4 , markerfacecolor='none', color='c' )
ax[4][1].plot(T2, abs(Uy2_err )/1e-16,  '<', label=r'$\epsilon_{U_y}$'  , ms=4 , markerfacecolor='none', color='b' )
ax[5][1].plot(T2, abs(Uz2_err )/1e-16,  'x', label=r'$\epsilon_{U_z}$'  , ms=4 , markerfacecolor='none', color='m' )

# Analytical error distribution
np.logspace( np.log10(min(T1)),np.log10(max(T1)), num=100 )
h1 = 2.5*T1 + np.sqrt(2.25*T1**2+1)
beta1 = np.sqrt( Ux1**2 + Uy1**2 + Uz1**2 ) / LorentzFactor1
AnalyticalError1  = (LorentzFactor1*h1)**2 * ( 1 + beta1**2 ) + (T1/LorentzFactor1)**2 - 2*h1*T1 - 1
AnalyticalError1 /= h1**2 + (T1/LorentzFactor1)**2 - 2*h1*T1 - 1
AnalyticalError1 *= np.finfo(np.float64).eps

ax[0][0].plot(T1, AnalyticalError1, '--',color='k',lw=3, label='Analytical') 
ax[1][0].plot(T1, AnalyticalError1, '--',color='k',lw=3) 
ax[2][0].plot(T1, AnalyticalError1, '--',color='k',lw=3) 
ax[3][0].plot(T1, AnalyticalError1, '--',color='k',lw=3) 
ax[4][0].plot(T1, AnalyticalError1, '--',color='k',lw=3) 
ax[5][0].plot(T1, AnalyticalError1, '--',color='k',lw=3) 

# Analytical error distribution
np.logspace( np.log10(min(T2)),np.log10(max(T2)), num=100 )
h2 = 2.5*T2 + np.sqrt(2.25*T2**2+1)
beta2 = np.sqrt( Ux2**2 + Uy2**2 + Uz2**2 ) / LorentzFactor2
AnalyticalError2  = (LorentzFactor2*h2)**2 * ( 1 + beta2**2 ) + (T2/LorentzFactor2)**2 - 2*h2*T2 - 1
AnalyticalError2 /= h2**2 + (T2/LorentzFactor2)**2 - 2*h2*T2 - 1
AnalyticalError2 *= np.finfo(np.float64).eps

ax[0][1].plot(T2, AnalyticalError2/1e-16, '--',color='k',lw=3, label='Analytical') 
ax[1][1].plot(T2, AnalyticalError2/1e-16, '--',color='k',lw=3) 
ax[2][1].plot(T2, AnalyticalError2/1e-16, '--',color='k',lw=3) 
ax[3][1].plot(T2, AnalyticalError2/1e-16, '--',color='k',lw=3) 
ax[4][1].plot(T2, AnalyticalError2/1e-16, '--',color='k',lw=3) 
ax[5][1].plot(T2, AnalyticalError2/1e-16, '--',color='k',lw=3) 


ax[5][0].set_xlabel('$k_{B}T/mc^2$'    , size='15')
ax[5][1].set_xlabel('$k_{B}T/mc^2$'    , size='15')

x_value = T1[np.where(np.isnan(Ux1_err))[0]]
ax[0][0].vlines(x_value,min(T1_err),max(T1_err),colors='b')



ax[0][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[1][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[2][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[3][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[4][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[5][1].set_ylabel('[$10^{-16}$]'    , size='15')




ax[0][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[1][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[2][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[3][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[4][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[5][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)

ax[0][0].set_yscale('log')
ax[1][0].set_yscale('log')
ax[2][0].set_yscale('log')
ax[3][0].set_yscale('log')
ax[4][0].set_yscale('log')
ax[5][0].set_yscale('log')


#ax[0][0].hlines(+np.finfo(np.float64).eps,min(T1),max(T1),colors='r', linestyle='--')
#ax[1][0].hlines(+np.finfo(np.float64).eps,min(T1),max(T1),colors='r', linestyle='--')
#ax[2][0].hlines(+np.finfo(np.float64).eps,min(T1),max(T1),colors='r', linestyle='--')
#ax[3][0].hlines(+np.finfo(np.float64).eps,min(T1),max(T1),colors='r', linestyle='--')
#ax[4][0].hlines(+np.finfo(np.float64).eps,min(T1),max(T1),colors='r', linestyle='--')
#ax[5][0].hlines(+np.finfo(np.float64).eps,min(T1),max(T1),colors='r', linestyle='--')

ax[0][0].set_xlim(min(T1),max(T1))
ax[1][0].set_xlim(min(T1),max(T1))
ax[2][0].set_xlim(min(T1),max(T1))
ax[3][0].set_xlim(min(T1),max(T1))
ax[4][0].set_xlim(min(T1),max(T1))
ax[5][0].set_xlim(min(T1),max(T1))

ax[0][1].set_xlim(min(T2),max(T2))
ax[1][1].set_xlim(min(T2),max(T2))
ax[2][1].set_xlim(min(T2),max(T2))
ax[3][1].set_xlim(min(T2),max(T2))
ax[4][1].set_xlim(min(T2),max(T2))
ax[5][1].set_xlim(min(T2),max(T2))


for axi in ax.flat:
  axi.tick_params(axis='both', labelsize='13', direction='in', bottom=True, top=True, right=True, left=True)
  axi.set_xscale('log')
  axi.minorticks_off()



ax[0][0].set_xticklabels([])
ax[1][0].set_xticklabels([])
ax[2][0].set_xticklabels([])
ax[3][0].set_xticklabels([])
ax[4][0].set_xticklabels([])

ax[0][1].set_xticklabels([])
ax[1][1].set_xticklabels([])
ax[2][1].set_xticklabels([])
ax[3][1].set_xticklabels([])
ax[4][1].set_xticklabels([])


# Set common labels
f.text(0.07, 0.5, 'relative error', ha='center', va='center', rotation='vertical', size=15)

# Set title
ax[0][0].set_title(r"$\gamma \sim 10^{%d}$" % np.log10(LorentzFactor1), size=15)
ax[0][1].set_title(r"$v/c \sim 10^{%d}$" % np.log10(beta2), size=15)



#plt.show()
plt.savefig('ErrorMap.pdf', format='pdf', bbox_inches='tight')
