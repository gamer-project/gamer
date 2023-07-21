import matplotlib
matplotlib.use('Agg')
import math
import h5py
import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end

def fit(t, a, b):
   return a*t+b

delta_t = 140.59 ## delta_t between snapshots (Myr)
t = []
for idx in range(idx_start, idx_end+1):
   t.append(idx*delta_t)
t = np.array(t)
mean_disp = np.zeros((4, len(t)))
i = 0
for idx in range(idx_start, idx_end+1):
   ParData = np.load('Vel_data_%06d.npz'%idx)
   mean_disp[0,i] = ParData['disp'][0]**2.
   mean_disp[1,i] = ParData['disp'][1]**2.
   mean_disp[2,i] = ParData['disp'][2]**2.
   mean_disp[3,i] = ParData['disp'][3]**2.
   i = i + 1
dt = delta_t*1e6*365*24*3600

heating = np.zeros(4)
const = np.zeros(4)
print('R(kpc)    heating rate (km^2/(s^2*Gyr))')
for i in range(4):
   xData = t
   yData = mean_disp[i,:]
   para_F, _ = curve_fit(fit, xData, yData)
   heating[i] = para_F[0]
   const[i] = para_F[1]
   print('%2d        %2.6f'%(i*2+4, heating[i]*1000 ))

plt.figure(dpi = 140)
plt.title('$\sigma_z^2$, bin size = 2 kpc')
plt.plot(t, mean_disp[0,:],'-o', ms=2, label = 'r = 4 kpc')
plt.plot(t, mean_disp[1,:],'-o', ms=2, label = 'r = 6 kpc')
plt.plot(t, mean_disp[2,:],'-o', ms=2, label = 'r = 8 kpc')
plt.plot(t, mean_disp[3,:],'-o', ms=2, label = 'r = 10 kpc')
for i in range(4):
   plt.plot(t, heating[i]*t+const[i], ms=2)
plt.grid(ls='--')
plt.legend(bbox_to_anchor=(1.04,1),loc='upper left', borderaxespad=0, shadow=True, prop={'size':8})   ## loc='best', 'upper left', 'upper right', 'lower left', 'lower right' 
plt.xlabel("t(Myr)")
plt.ylabel("$\sigma_z^2\, (km/s)$")
plt.savefig("sigma_z_sqr.png", dpi = 140, bbox_inches="tight")
plt.close()


