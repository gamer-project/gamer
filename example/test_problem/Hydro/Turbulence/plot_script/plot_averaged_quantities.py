import argparse
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np
import gc

# load the command-line parameters
parser = argparse.ArgumentParser( description='Compute averaged quantities' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start    = args.idx_start
idx_end      = args.idx_end
didx         = args.didx

dpi          = 150

yt.enable_parallelism()


# load the dataset
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
my_storage = {}

# main loop
for sto, ds in ts.piter(storage=my_storage):

   N  = np.int64( ds.parameters["NX0"][0] )
   dd = ds.covering_grid( level=0, left_edge=[0, 0, 0], dims=ds.domain_dimensions )
   dens = dd["Dens"].d

#  remove center of mass motion
   vx = dd["MomX"].d / dens - np.mean( dd['MomX'] ).d / np.mean( dens )
   vy = dd["MomY"].d / dens - np.mean( dd['MomY'] ).d / np.mean( dens )
   vz = dd["MomZ"].d / dens - np.mean( dd['MomZ'] ).d / np.mean( dens )

#  magnetic field
   Bx = dd['magnetic_field_x'].d
   By = dd['magnetic_field_y'].d
   Bz = dd['magnetic_field_z'].d

#  mean magnetic_field
   mean_Bx = np.mean( Bx )
   mean_By = np.mean( By )
   mean_Bz = np.mean( Bz )

   kx = (2.0 * np.pi * np.fft.fftfreq (N, d=1.0/N))[:, None, None]
   ky = (2.0 * np.pi * np.fft.fftfreq (N, d=1.0/N))[None, :, None]
   kz = (2.0 * np.pi * np.fft.rfftfreq(N, d=1.0/N))[None, None, :]
   kk = kx**2 + ky**2 + kz**2
   kk[0,0,0] = 1.0
   weight = np.ones_like(kk, dtype=float)
   weight[:, :, 1:-1] = 2.0

#  vorticity w = curl(v) -> wk = (ik x v)
   vxk  = np.fft.rfftn( vx )
   vyk  = np.fft.rfftn( vy )
   vzk  = np.fft.rfftn( vz )

   wxk = 1j * (ky * vzk - kz * vyk)
   wyk = 1j * (kz * vxk - kx * vzk)
   wzk = 1j * (kx * vyk - ky * vxk)

#  averaged kinetic helicity = sum(w dot v)/N^3 = sum(wk dot vk)/N^6
   hkin = np.sum( weight * (np.conj(wxk)*vxk + np.conj(wyk)*vyk + np.conj(wzk)*vzk) ).real / N**6

   del vxk, vyk, vzk, wxk, wyk, wzk
   gc.collect()

   Bxk  = np.fft.rfftn( Bx )
   Byk  = np.fft.rfftn( By )
   Bzk  = np.fft.rfftn( Bz )

#  vector potential (without DC): Laplace(A) = -curl(B) (Coulomb gauge) -> Ak = (ik x B) / k^2
   Axk = 1j * (ky * Bzk - kz * Byk) / kk
   Ayk = 1j * (kz * Bxk - kx * Bzk) / kk
   Azk = 1j * (kx * Byk - ky * Bxk) / kk

#  current J = curl(B) -> Jk = (ik x B)
   Jxk = 1j * (ky * Bzk - kz * Byk)
   Jyk = 1j * (kz * Bxk - kx * Bzk)
   Jzk = 1j * (kx * Byk - ky * Bxk)

#  averaged magnetic helicity (fluctuated) = sum(A dot B)/N^3 = sum(Ak dot Bk)/N^6
   hmag = np.sum( weight * (np.conj(Axk)*Bxk + np.conj(Ayk)*Byk + np.conj(Azk)*Bzk) ).real / N**6

#  averaged current helicity = sum(J dot B)/N^3 = sum(Jk dot Bk)/N^6
   hcur = np.sum( weight * (np.conj(Jxk)*Bxk + np.conj(Jyk)*Byk + np.conj(Jzk)*Bzk) ).real / N**6

   del Bxk, Byk, Bzk, Axk, Ayk, Azk, Jxk, Jyk, Jzk
   gc.collect()

   emag   = np.mean( dd['magnetic_energy_density'].in_units('code_mass*code_velocity**2/code_length**3').d )
   ek_tot = np.mean( dd['kinetic_energy_density' ].in_units('code_mass*code_velocity**2/code_length**3').d )
   ekin   = np.mean( 0.5*dens*( vx*vx + vy*vy + vz*vz ) )
   mach   = np.mean( dd['mach_number'].d )

   time = ds.current_time
   sto.result = {
        "time"   : time,
        "emag"   : emag,
        "ekin"   : ekin,
        "ek_tot" : ek_tot,
        "mach"   : mach,
        "hkin"   : hkin,
        "hmag"   : hmag,
        "hcur"   : hcur,
        "meanBx" : mean_Bx,
        "meanBy" : mean_By,
        "meanBz" : mean_Bz,
   }

# plot
if yt.is_root():
   fields = [
       "time",
       "emag",
       "ekin",
       "ek_tot",
       "mach",
       "hkin",
       "hmag",
       "hcur",
       "meanBx",
       "meanBy",
       "meanBz",
   ]

   data = { field: np.array([val[field] for dx, val in sorted(my_storage.items())]) for field in fields }

   time   = data["time"  ]
   emag   = data["emag"  ]
   ekin   = data["ekin"  ]
   ek_tot = data["ek_tot"]
   mach   = data["mach"  ]
   hkin   = data["hkin"  ]
   hmag   = data["hmag"  ]
   hcur   = data["hcur"  ]
   meanBx = data["meanBx"]
   meanBy = data["meanBy"]
   meanBz = data["meanBz"]

   ratio = np.divide( emag, ekin,
           out=np.zeros_like(emag, dtype=float),
           where=ekin != 0 )
   # Plot
   f, ax = plt.subplots(2, 1, figsize=(6, 8))
   f.subplots_adjust(wspace=0.4)
   ax[0].plot( time[1:], ratio[1:], label = r"$E_{\rm mag}/E_{\rm kin}$" )
   ax[0].set_yscale('log')
   ax[0].set_xlim(0, 100)
   ax[0].set_ylabel( r'$E_{\rm mag}/E_{\rm kin}$', fontsize='large' )

   ax[1].plot( time, mach, label = "Mach number")
   ax[1].set_xlim(0, 100)
   ax[1].set_ylim(0, 0.4)
   ax[1].set_ylabel('Mach number', fontsize='large')
   ax[1].set_xlabel('time',        fontsize='large')

   plt.savefig("fig_EnergyEvolve.png", bbox_inches='tight', pad_inches=0.05, dpi=dpi)
   plt.close()

   # save figure
   np.savetxt( 'AveragedQuantities', np.column_stack( (time, emag, ekin, ek_tot, mach, hkin, hmag, hcur, meanBx, meanBy, meanBz) ),
               fmt='  %16.8e',
               header='               t               Emag               Ekin             Ek_tot               Mach               Hkin               Hmag               Hcur            mean_Bx            mean_By            mean_Bz' )


