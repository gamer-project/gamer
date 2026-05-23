import numpy as np
import matplotlib.pyplot as plt
import yt
import argparse
import sys

#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Get power spectrum' )

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

# plot GAMER generated power spectrum as reference
GAMER_PS = True
print('GAMER_PS='+str(GAMER_PS)+'\n')

yt.enable_parallelism()
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
   idx     = ds.parameters["DumpID"]
   time    = ds.parameters["Time"][0]
   N       = ds.parameters["NX0"][0]
   BoxSize = ds.parameters["BoxSize"][0]
   dh      = ds.parameters["CellSize"][0]
   dd      = ds.covering_grid(level=0, left_edge=[0, 0, 0], dims=ds.domain_dimensions)
   rho     = dd["Dens"].d
   rhok    = np.fft.rfftn( rho )

   kx = 2*np.pi*np.fft.fftfreq( N, dh )
   ky = 2*np.pi*np.fft.fftfreq( N, dh )
   kz = 2*np.pi*np.fft.rfftfreq( N, dh )

   Pk3d = abs(rhok)**2
   Pk_total = np.zeros(N//2+1)
   count = np.zeros(N//2+1)

   for i in range( N ):
      for j in range( N ):
         for k in range( N//2+1 ):
            l = round((kx[i]**2+ky[j]**2+kz[k]**2)**0.5*BoxSize/(2.0*np.pi))
            if (l < N//2+1):
               Pk_total[l] = Pk_total[l] + Pk3d[i,j,k]
               count[l] = count[l] + 1
   Pk1d = np.divide(Pk_total, count, out=np.zeros_like(Pk_total), where=count!=0)
   k3Pk = Pk1d*kz**3
   d = 2*np.pi/kz[np.argmax(k3Pk)]
   print('estimated granule diameter = %13.7e, time = %13.7e\n'%(d, time))

   if(GAMER_PS):
      gamer_ps = np.loadtxt('../PowerSpec_%06d'%idx, skiprows=1, dtype=float)
      gamer_k    = gamer_ps[:,0]
      gamer_Pk1d = gamer_ps[:,1]
      gamer_k3Pk = gamer_k**3*gamer_Pk1d

   plt.figure()
   plt.title("Dimensionless Power Spectrum")
   plt.plot(kz[1:], k3Pk[1:]/max(k3Pk), label = 'numpy')
   if(GAMER_PS):
      plt.plot(gamer_k, gamer_k3Pk/max(gamer_k3Pk), label = 'gamer')
   plt.xlabel('$k$')
   plt.ylabel(r'$k^3P(k)$')
   plt.yscale('log')
   plt.legend(loc='lower right')
   plt.savefig('fig_powerspectrum_dimensionless_%06d.png'%idx, dpi = 150, bbox_inches="tight")
   plt.close()


