import argparse
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Particle density profile' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
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
prefix       = args.prefix

r_sphere     = (100.0, 'kpc')
dpi          = 150
nbin         = 32


# add the total density field, including the fluid density and the particle density
def _TotDens(field, data):
    return data['Dens']+data['ParDens']

yt.add_field( ('gamer', 'TotDens'),
              function=_TotDens, units='code_mass/code_length**3',
              sampling_type='cell' )

# add the total mass field, including the fluid mass and the particle mass
def _TotCellMass(field, data):
    return data['TotDens']*data['cell_volume']

yt.add_field( ('gamer', 'TotCellMass'),
              function=_TotCellMass, units='code_mass',
              sampling_type='cell' )


yt.enable_parallelism()
ts = yt.DatasetSeries([ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ])

for ds in ts.piter():
   for field_mass in [('gamer', 'TotCellMass')]:

      # set the center as the position of the maximum density
      center_pos  = ds.all_data().quantities.max_location(('gamer', 'TotDens'))[1:]
      sp          = ds.sphere( center_pos, r_sphere )

      # create the arrays of profile
      prof_mass = yt.create_profile( sp, 'radius', fields=field_mass,
                                     weight_field=None, n_bins=nbin, accumulation=True,
                                     units={'radius':'code_length', field_mass:'code_mass'} )

      # save the profile to text file
      np.savetxt( '%s_EnclosedMass_profile'%(ds),
                  np.column_stack( (prof_mass.x.in_units('code_length').d, prof_mass[field_mass].in_units('code_mass').d) ),
                  fmt='          %9.8e',
                  header='       r (code_length)         mass (code_mass)' )

      # decide the units for plotting
      #UNIT_L_PLOT = 'code_length'
      #UNIT_M_PLOT = 'code_mass'
      UNIT_L_PLOT = 'kpc'
      UNIT_M_PLOT = 'Msun'

      # create the figure
      fig = plt.figure()
      ax  = fig.add_subplot(111)

      # plot the profiles
      ax.plot( prof_mass.x.in_units(UNIT_L_PLOT).d, prof_mass[field_mass].in_units(UNIT_M_PLOT).d,  label='%s'%ds )

      # set the limits and scales
      ax.set_xlim( 1.0e-1, 1.0e+2 )
      ax.set_ylim( 1.0e+5, 1.0e10 )
      ax.set_xscale('log')
      ax.set_yscale('log')

      # set the labels
      ax.set_xlabel( r'$r$'+' (%s)'%UNIT_L_PLOT )
      ax.set_ylabel( r'$M$'+' (%s)'%UNIT_M_PLOT )
      ax.legend()
      fig.suptitle( '$t$ = %7.6e Gyr'%ds.current_time.in_units('Gyr') )
      ax.annotate( 'Total Mass = %7.6e Msun'%ds.sphere( center_pos, r_sphere ).quantities.total_quantity([field_mass]).in_units('Msun'),
                   xy=(0.5,0.05), xycoords='axes fraction' )

      # set the grid and ticks
      ax.grid()
      ax.xaxis.set_ticks_position('both')
      ax.yaxis.set_ticks_position('both')
      ax.tick_params( which='both',direction='in' )

      # save the figure
      plt.tight_layout( pad=0.1, h_pad=0.1, w_pad=0.1 )
      fig.savefig( '%s_EnclosedMass_profile.png'%(ds), dpi=dpi )
      fig.clear()
