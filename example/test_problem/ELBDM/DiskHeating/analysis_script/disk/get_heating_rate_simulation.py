import matplotlib
matplotlib.use('Agg')
import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.pyplot import cm
from matplotlib.patches import Rectangle


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Compute the time-averaged simulating heating rate' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )


#-----------------------------------------------------------------------------------------------------------------
# predefined functions
def fit(t, a, b):
   return a*t+b


#-----------------------------------------------------------------------------------------------------------------
# compute disk heating rates from simulation data
delta_t   = 140.59           # delta_t between snapshots (Myr)
dt        = delta_t*1e6*365*24*3600
t         = []
for idx in range(idx_start, idx_end+1):
   t.append(idx*delta_t)
t         = np.array(t)
mean_disp = np.zeros((4, len(t)))
i         = 0
for idx in range(idx_start, idx_end+1):
   ParData = np.load('Vel_data_%06d.npz'%idx)
   mean_disp[0,i] = ParData['disp'][0]**2.
   mean_disp[1,i] = ParData['disp'][1]**2.
   mean_disp[2,i] = ParData['disp'][2]**2.
   mean_disp[3,i] = ParData['disp'][3]**2.
   i = i + 1

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

# plot result
lines=[]
color = [cm.Blues(np.linspace(0.4, 0.7, 2)),
         cm.Oranges(np.linspace(0.3, 0.6, 2)),
         cm.Greens(np.linspace(0.4, 0.7, 2)),
         cm.Reds(np.linspace(0.5, 0.8, 2))]

fig = plt.figure()
axs = fig.add_subplot(111)
plt.title('$\sigma_z^2$, bin size = 2 kpc')
for i in range(4):
  l1, = axs.plot(t/1000., mean_disp[i,:],'-o', ms=2, c=color[i][1])
  l2, = axs.plot(t/1000., heating[i]*t+const[i],'--', dashes=(5,1), ms=2, c=color[i][0])
  lines.append([l1, l2])

# create legend
extra = Rectangle((0, 0), 1, 1.0, fc="w", fill=False, edgecolor='none', linewidth=0)

legend_handle = [extra,    extra   ,    extra   ,    extra   ,    extra,
                 extra, lines[0][0], lines[1][0], lines[2][0], lines[3][0],
                 extra, lines[0][1], lines[1][1], lines[2][1], lines[3][1]]

length = len(legend_handle)
legend_handle = np.array(legend_handle).reshape(int(length/5), 5).T.reshape(length)

label_row1 = [''          , '$R$ = 4 kpc', '$R$ = 6 kpc', '$R$ = 8 kpc', '$R$ = 10 kpc']
label_row2 = ['simulation', ''           , ''           , ''           , ''            ]
label_row3 = ['linear fit', ''           , ''           , ''           , ''            ]
legend_labels = np.concatenate([label_row1, label_row2, label_row3]).reshape(int(length/5), 5).T.reshape(length)

axs.legend(legend_handle, legend_labels, loc='upper left',fontsize=7, ncol = 5, handletextpad = -2.5, handlelength=2.5, handleheight=1.5)

plt.ylim(0, 1.1*max(mean_disp[0]))
plt.xlabel(r"$t_{\rm rel}$ (Gyr)")
plt.ylabel("$\sigma_z^2$ (km$^2$/s$^2$)")
plt.savefig("sigma_z_sqr.png", dpi = 140, bbox_inches="tight")
plt.close()
