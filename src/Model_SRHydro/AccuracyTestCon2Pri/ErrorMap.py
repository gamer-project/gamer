import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import os

pwd = os.getcwd()
pwd = pwd.split('/')

f, ax = plt.subplots( 8, 2, sharex=False, sharey=False )
f.subplots_adjust( hspace=0.05, wspace=0.15 )
f.set_size_inches( 12.0, 10.0 )

def MachNumber_high(T):
  h = 2.5*T + np.sqrt( 2.25 * T**2 + 1.0 ) 
  a = T *( 5.0*h - 8.0*T )
  b = 3.0 * h * (h - T)
  Cs_sq = a / b 
  Cs = np.sqrt( Cs_sq )
  Us = Cs / np.sqrt( 1.0 - Cs_sq )
  return 1e-4/Us

def MachNumber_low(T):
  h = 2.5*T + np.sqrt( 2.25 * T**2 + 1.0 ) 
  a = T *( 5.0*h - 8.0*T )
  b = 3.0 * h * (h - T)
  Cs_sq = a / b 
  Cs = np.sqrt( Cs_sq )
  Us = Cs / np.sqrt( 1.0 - Cs_sq )
  return 1e-4/Us
 

def _HTilde (T):
  HTilde  = 2.5*T + 2.25*T**2
  HTilde /= np.sqrt(2.25*T**2+1) + 1
  return HTilde

X1, Y1 = [], []
X2, Y2 = [], []
data = []

data1 = np.loadtxt('ErrorMapURlimit.dat')
data2 = np.loadtxt('ErrorMapNRlimit.dat')

rho1_err       = data1[:,0]  
Ux1_err        = data1[:,1]  
Uy1_err        = data1[:,2]  
Uz1_err        = data1[:,3]  
P1_err         = data1[:,4]  
T1_err         = data1[:,5]  
HT1_err        = data1[:,6]
H1_err         = data1[:,7]

Ux1            = data1[0,10]
Uy1            = data1[0,11]
Uz1            = data1[0,12]
LorentzFactor1 = data1[0,14]
T1             = data1[:,15]  
HT1            = data1[:,16]


rho2_err       = data2[:,0]  
Ux2_err        = data2[:,1]  
Uy2_err        = data2[:,2]  
Uz2_err        = data2[:,3]  
P2_err         = data2[:,4]  
T2_err         = data2[:,5]  
HT2_err        = data2[:,6]
H2_err         = data2[:,7]

Ux2            = data2[0,10]
Uy2            = data2[0,11]
Uz2            = data2[0,12]
LorentzFactor2 = data2[0,14]
T2             = data2[:,15]  
HT2            = data2[:,16]

ax[0][0].plot(T1, abs(HT1_err ),  'x', label=r'$\tilde{h}$'  , ms=4 , markerfacecolor='none', color='red'    )
ax[1][0].plot(T1, abs(T1_err  ),  '+', label=r'$T$'          , ms=4 , markerfacecolor='none', color='salmon' )
ax[2][0].plot(T1, abs(P1_err  ),  'D', label=r'$P$'          , ms=4 , markerfacecolor='none', color='pink'   )
ax[3][0].plot(T1, abs(H1_err  ),  '1', label=r'$h$'          , ms=4 , markerfacecolor='none', color='orange' )
ax[4][0].plot(T1, abs(rho1_err),  'o', label=r'$\rho$'       , ms=4 , markerfacecolor='none', color='gold'   )
ax[5][0].plot(T1, abs(Ux1_err ),  '^', label=r'$U_x$'        , ms=4 , markerfacecolor='none', color='yellowgreen'   )
ax[6][0].plot(T1, abs(Uy1_err ),  '<', label=r'$U_y$'        , ms=4 , markerfacecolor='none', color='orchid' )
ax[7][0].plot(T1, abs(Uz1_err ),  '2', label=r'$U_z$'        , ms=4 , markerfacecolor='none', color='fuchsia')

ax[0][1].plot(T2, abs(HT2_err )/1e-16,  'x', label=r'$\tilde{h}$', ms=4 , markerfacecolor='none', color='red'    )
ax[1][1].plot(T2, abs(T2_err  )/1e-16,  '+', label=r'$T$'        , ms=4 , markerfacecolor='none', color='salmon' )
ax[2][1].plot(T2, abs(P2_err  )/1e-16,  'D', label=r'$P$'        , ms=4 , markerfacecolor='none', color='pink'   )
ax[3][1].plot(T2, abs(H2_err  )/1e-16,  'x', label=r'$h$'        , ms=4 , markerfacecolor='none', color='orange' )
ax[4][1].plot(T2, abs(rho2_err)/1e-16,  'o', label=r'$\rho$'     , ms=4 , markerfacecolor='none', color='gold'   )
ax[5][1].plot(T2, abs(Ux2_err )/1e-16,  '^', label=r'$U_x$'      , ms=4 , markerfacecolor='none', color='yellowgreen'   )
ax[6][1].plot(T2, abs(Uy2_err )/1e-16,  '<', label=r'$U_y$'      , ms=4 , markerfacecolor='none', color='orchid' )
ax[7][1].plot(T2, abs(Uz2_err )/1e-16,  'x', label=r'$U_z$'      , ms=4 , markerfacecolor='none', color='fuchsia')

# Analytical error distribution
np.logspace( np.log10(min(T1)),np.log10(max(T1)), num=100 )
h1 = 2.5*T1 + np.sqrt(2.25*T1**2+1)
beta1 = np.sqrt( Ux1**2 + Uy1**2 + Uz1**2 ) / LorentzFactor1
AnalyticalError1  = (LorentzFactor1*h1)**2 * ( 1 + beta1**2 ) + (T1/LorentzFactor1)**2 - 2*h1*T1 - 1
AnalyticalError1 /= h1**2 + (T1/LorentzFactor1)**2 - 2*h1*T1 - 1
AnalyticalError1 *= np.finfo(np.float64).eps


AnalyticalError1_U = ( HT1 - _HTilde(T1) + HT1*AnalyticalError1 ) / ( 1 + HT1 + HT1*AnalyticalError1 )

ax[0][0].plot(T1, AnalyticalError1,            ':', color='black',lw=3) 
ax[1][0].plot(T1, AnalyticalError1,            ':', color='black',lw=3) 
ax[2][0].plot(T1, AnalyticalError1,            ':', color='black',lw=3) 
ax[3][0].plot(T1, AnalyticalError1*HT1/(1+HT1),'--',color='black',lw=3) 
ax[4][0].plot(T1, AnalyticalError1*HT1/(1+HT1),'--',color='black',lw=3) 
ax[5][0].plot(T1, AnalyticalError1*HT1/(1+HT1),'--',color='black',lw=3) 
ax[6][0].plot(T1, AnalyticalError1*HT1/(1+HT1),'--',color='black',lw=3) 
ax[7][0].plot(T1, AnalyticalError1*HT1/(1+HT1),'--',color='black',lw=3) 

# Analytical error distribution
np.logspace( np.log10(min(T2)),np.log10(max(T2)), num=100 )
h2 = 2.5*T2 + np.sqrt(2.25*T2**2+1)
beta2 = np.sqrt( Ux2**2 + Uy2**2 + Uz2**2 ) / LorentzFactor2
AnalyticalError2  = (LorentzFactor2*h2)**2 * ( 1 + beta2**2 ) + (T2/LorentzFactor2)**2 - 2*h2*T2 - 1
AnalyticalError2 /= h2**2 + (T2/LorentzFactor2)**2 - 2*h2*T2 - 1
AnalyticalError2 *= np.finfo(np.float64).eps

ax[0][1].plot(T2, AnalyticalError2/1e-16,              '-',color='black',lw=3) 
ax[1][1].plot(T2, AnalyticalError2/1e-16,              '-',color='black',lw=3) 
ax[2][1].plot(T2, AnalyticalError2/1e-16,              '-',color='black',lw=3) 
ax[3][1].plot(T2, AnalyticalError2*HT2/(1+HT2)/1e-16, '-.',color='black',lw=3) 
ax[4][1].plot(T2, AnalyticalError2*HT2/(1+HT2)/1e-16, '-.',color='black',lw=3) 
ax[5][1].plot(T2, AnalyticalError2*HT2/(1+HT2)/1e-16, '-.',color='black',lw=2) 
ax[6][1].plot(T2, AnalyticalError2*HT2/(1+HT2)/1e-16, '-.',color='black',lw=3) 
ax[7][1].plot(T2, AnalyticalError2*HT2/(1+HT2)/1e-16, '-.',color='black',lw=3) 


ax[7][0].set_xlabel('$k_{B}T/mc^2$'    , size='15')
ax[7][1].set_xlabel('$k_{B}T/mc^2$'    , size='15')

x_value = T1[np.where(np.isnan(Ux1_err))[0]]
ax[0][0].vlines(x_value,min(T1_err),max(T1_err),colors='b')



ax[0][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[1][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[2][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[3][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[4][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[5][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[6][1].set_ylabel('[$10^{-16}$]'    , size='15')
ax[7][1].set_ylabel('[$10^{-16}$]'    , size='15')




ax[0][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[1][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[2][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[3][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[4][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[5][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[6][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)
ax[7][0].legend(loc='upper right', markerscale=1, handletextpad=0.1,columnspacing=2,ncol=3, fontsize=15, borderaxespad=0.1)

ax[0][0].set_yscale('log')
ax[1][0].set_yscale('log')
ax[2][0].set_yscale('log')
ax[3][0].set_yscale('log')
ax[4][0].set_yscale('log')
ax[5][0].set_yscale('log')
ax[6][0].set_yscale('log')
ax[7][0].set_yscale('log')

ax[0][0].set_ylim(1e-6,5e0)
ax[1][0].set_ylim(1e-6,5e0)
ax[2][0].set_ylim(1e-6,5e0)
ax[3][0].set_ylim(1e-6,5e0)
ax[4][0].set_ylim(1e-6,5e0)
ax[5][0].set_ylim(1e-6,5e0)
ax[6][0].set_ylim(1e-6,5e0)
ax[7][0].set_ylim(1e-6,5e0)

ax[0][1].set_ylim(-1.5,8.0)
ax[1][1].set_ylim(-1.5,8.0)
ax[2][1].set_ylim(-1.5,8.0)
ax[3][1].set_ylim(-1.5,8.0)
ax[4][1].set_ylim(-1.5,8.0)
ax[5][1].set_ylim(-1.5,8.0)
ax[6][1].set_ylim(-1.5,8.0)
ax[7][1].set_ylim(-1.5,8.0)


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
ax[6][0].set_xlim(min(T1),max(T1))
ax[7][0].set_xlim(min(T1),max(T1))

ax[0][1].set_xlim(min(T2),max(T2))
ax[1][1].set_xlim(min(T2),max(T2))
ax[2][1].set_xlim(min(T2),max(T2))
ax[3][1].set_xlim(min(T2),max(T2))
ax[4][1].set_xlim(min(T2),max(T2))
ax[5][1].set_xlim(min(T2),max(T2))
ax[6][1].set_xlim(min(T2),max(T2))
ax[7][1].set_xlim(min(T2),max(T2))


for axi in ax.flat:
  axi.tick_params(axis='both', labelsize='13', direction='in', bottom=True, top=True, right=True, left=True)
  axi.set_xscale('log')
  axi.minorticks_off()


ax[0][0].set_xticklabels([])
ax[1][0].set_xticklabels([])
ax[2][0].set_xticklabels([])
ax[3][0].set_xticklabels([])
ax[4][0].set_xticklabels([])
ax[5][0].set_xticklabels([])
ax[6][0].set_xticklabels([])

ax[0][1].set_xticklabels([])
ax[1][1].set_xticklabels([])
ax[2][1].set_xticklabels([])
ax[3][1].set_xticklabels([])
ax[4][1].set_xticklabels([])
ax[5][1].set_xticklabels([])
ax[6][1].set_xticklabels([])


# Set common labels
f.text(0.07, 0.5, 'relative error', ha='center', va='center', rotation='vertical', size=15)

# Set title
ax[0][0].set_title(r"$\gamma = 10^{%d}$" % np.log10(LorentzFactor1), size=15)
ax[0][1].set_title(r"$\beta = 10^{%d}$" % np.log10(beta2), size=15)



plt.show()
#plt.savefig('ErrorMap.pdf', format='pdf', bbox_inches='tight')
