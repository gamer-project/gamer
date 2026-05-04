import matplotlib
matplotlib.use('Agg')
import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
output_mode = int(3)         # [1] angle-averaged rotation curve; [2] shell-averaged halo density profile; [3] both
figure_dpi  = 140


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='An example script to plot disk/halo data from existing .npy files' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )


#-------------------------------------------------------------------------------------------------------------------------
# analyze and plot input simulation data
# [1] angle-averaged rotation curve
if (output_mode == 1 or output_mode == 3):
   plt.figure(dpi = figure_dpi)
   for idx in range(idx_start, idx_end+1, didx):

      data = np.load('disk/Data_Disk_%06d.npy'%idx)
      r = data[0]
      r_in_kpc = r/3.08568e+21
      v_c = data[1]/1.0e5       # rotaion speed
      sigma_r = data[2]/1.0e5   # sigma_r
      sigma_z = data[3]/1.0e5   # sigma_z
      Sigma = data[4]           # surface density
      height = data[7]          # scale height
      plt.plot(r_in_kpc, v_c, color = cm.Blues(0.3+0.6*(idx-idx_start)/(idx_end+1-idx_start)), label ='id=%d'%idx)

   plt.xlim((0, 15))
   plt.ylim((0, 100))
   plt.grid(ls='--')
   plt.legend(loc='lower right', shadow=True, prop={'size':8})   ## loc='best', 'upper left', 'upper right', 'lower left', 'lower right'
   plt.xlabel(r"$R$ (kpc)")
   plt.ylabel(r"$v_{\rm cir}$ (km/s)")
   plt.savefig("Rotation_Curve.png", dpi = figure_dpi)
   plt.close()

# [2] plot halo density profile
if (output_mode == 2 or output_mode == 3):
   plt.figure(dpi = figure_dpi)
   for idx in range(idx_start, idx_end+1, didx):

      Dens =  np.load('halo/Halo_Dens_Data_%06d.npy'%idx)
      Dens =  np.delete(Dens,  np.nonzero(Dens[1]==0)[0] , axis = 1)
      r_d = Dens[0]
      dens_cgs = Dens[1]
      plt.plot(r_d, dens_cgs, color = cm.Blues(0.3+0.6*(idx-idx_start)/(idx_end+1-idx_start)), label ='id=%d'%idx)

   plt.xlim((1e-1, 1e2))
   plt.grid(ls='--')
   plt.legend(loc='lower left', shadow=True, prop={'size':8})   ## loc='best', 'upper left', 'upper right', 'lower left', 'lower right'
   plt.xlabel(r"$r$ (kpc)")
   plt.ylabel(r"$\rho_{\rm h}$ (g/cm$^3$)")
   plt.xscale('log')
   plt.yscale('log')
   plt.savefig("Halo_Density_Profile.png", dpi = figure_dpi)
   plt.close()
