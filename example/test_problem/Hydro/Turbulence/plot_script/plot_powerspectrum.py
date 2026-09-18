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

yt.enable_parallelism()
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
   idx     = ds.parameters["DumpID"]
   time    = ds.parameters["Time"][0]
   BoxSize = ds.parameters["BoxSize"][0]
   dh      = ds.parameters["CellSize"][0]
   kmax    = ds.parameters["Src_Turb_Kmax"]
   N       = np.int64( ds.parameters["NX0"][0] )
   dd      = ds.covering_grid(level=0, left_edge=[0, 0, 0], dims=ds.domain_dimensions)
   dens    = dd["Dens"].d
   wx      = dd["MomX"].d/(dens)**0.5
   wy      = dd["MomY"].d/(dens)**0.5
   wz      = dd["MomZ"].d/(dens)**0.5
   bx      = dd['CCMagX'].d
   by      = dd['CCMagY'].d
   bz      = dd['CCMagZ'].d

   wxk  = np.fft.rfftn( wx )
   wyk  = np.fft.rfftn( wy )
   wzk  = np.fft.rfftn( wz )
   bxk  = np.fft.rfftn( bx )
   byk  = np.fft.rfftn( by )
   bzk  = np.fft.rfftn( bz )

   Pk3d_kin = 0.5*( np.abs(wxk)**2 + np.abs(wyk)**2 + np.abs(wzk)**2 )
   Pk3d_mag = 0.5*( np.abs(bxk)**2 + np.abs(byk)**2 + np.abs(bzk)**2 )

   Pk_total_kin = np.zeros(N//2+1)
   Pk_total_mag = np.zeros(N//2+1)

   ix, iy, iz = np.meshgrid( np.fft.fftfreq(N)*N, np.fft.fftfreq(N)*N, np.fft.rfftfreq(N)*N, indexing="ij" )

   l = np.rint(np.sqrt(ix**2 + iy**2 + iz**2)).astype(int)

   weights = np.ones_like(l, dtype=float)
   weights[:, :, 1:-1] = 2.0

#  compute total energy of each k shell
   Pk_total_kin = np.bincount( l.ravel(), weights=(weights * Pk3d_kin).ravel(), minlength=N//2+1 )[:N//2+1]
   Pk_total_mag = np.bincount( l.ravel(), weights=(weights * Pk3d_mag).ravel(), minlength=N//2+1 )[:N//2+1]

#  energy density power spectrum normalization
#  Parseval's theorem: Sum(v^2)/N^3 = Sum(vk^2)/N^6
   Pk_total_kin /= N**6
   Pk_total_mag /= N**6

#  check
   print("\n  t = %f :\n"%time
         + "  average kinetic  energy density in x space = %13.7e\n"%(np.sum(wx*wx+wy*wy+wz*wz)*0.5/N**3)
         + "  average kinetic  energy density in k space = %13.7e\n"%(np.sum(Pk_total_kin))
         + "  average magnetic energy density in x space = %13.7e\n"%(np.sum(bx*bx+by*by+bz*bz)*0.5/N**3)
         + "  average magnetic energy density in k space = %13.7e\n"%(np.sum(Pk_total_mag))
        )

   k  = np.arange(N//2+1) * (2*np.pi/BoxSize)
   dk = k[1] - k[0]

#  divide by dk so Sum(Pk*dk) is the total energy
   Pk_total_kin /= dk
   Pk_total_mag /= dk

   plt.figure()
   plt.title("Power Spectrum")
   plt.plot(k[1:], Pk_total_kin[1:], label = r'$E_{\rm kin}(k)$')
   plt.plot(k[1:], Pk_total_mag[1:], label = r'$E_{\rm mag}(k)$')
   plt.xlabel('$k$')
   plt.ylabel(r'$E(k)$')
   plt.yscale('log')
   plt.xscale('log')
   plt.axvline(kmax*2*np.pi/BoxSize, color = '0.8', ls = '--')
   plt.legend(loc='upper right')
   plt.savefig('fig_powerspectrum_%06d.png'%idx, dpi = 150, bbox_inches="tight")
   plt.close()

   np.savetxt( 'EnergyPowerSpec_%06d'%idx, np.column_stack( (k, Pk_total_mag, Pk_total_kin)),
               fmt='  %16.8e',
               header='               k               EMag               Ekin' )

