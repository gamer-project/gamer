import argparse
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Profiles' )

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

r_sphere     = (1.5, 'code_length')
dpi          = 150
nbin         = 64


yt.enable_parallelism()
ts = yt.DatasetSeries([ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ])

for ds in ts.piter():

   fields_list = ['particle_density_on_grid']

   for field in fields_list:

      # set the center as the position of the gravitational potential minimum
      center_pos  = ds.all_data().quantities.max_location(('gas','particle_density_on_grid'))[1:]
      sp          = ds.sphere( center_pos, r_sphere )

      # create the arrays of profile
      prof_dens = yt.create_profile( sp, 'radius', fields=field,
                                     weight_field='cell_volume', n_bins=nbin,
                                     units={'radius':'code_length', field:'code_density'} )

      # save the profile to text file
      np.savetxt( '%s_%s_profile'%(ds,field),
                  np.column_stack( (prof_dens.x.in_units('code_length').d, prof_dens[field].in_units('code_density').d) ),
                  fmt='          %9.8e',
                  header='       r (code_length)   density (code_density)' )

      # load the reference profiles
      Ini_DensProf_r, Ini_DensProf_dens = np.loadtxt( '%s_%s_profile'%('Data_000000',field), skiprows=1, unpack=True )

      # assign the units
      Ini_DensProf_r    = ds.arr( Ini_DensProf_r,    'code_length'  )
      Ini_DensProf_dens = ds.arr( Ini_DensProf_dens, 'code_density' )

      # decide the units for plotting
      UNIT_L_PLOT = 'code_length'
      UNIT_D_PLOT = 'code_density'

      # create the figure
      fig = plt.figure()
      ax  = fig.add_subplot(111)

      # plot the profiles
      ax.plot( Ini_DensProf_r.in_units(UNIT_L_PLOT).d, Ini_DensProf_dens.in_units(UNIT_D_PLOT).d, label='Data_000000' )
      ax.plot( prof_dens.x.in_units(UNIT_L_PLOT).d,    prof_dens[field].in_units(UNIT_D_PLOT).d,  label='%s'%ds       )

      # set the limits and scales
      ax.set_xlim( 0.5*np.min(Ini_DensProf_r.in_units(UNIT_L_PLOT).d),     2.0*np.max(Ini_DensProf_r.in_units(UNIT_L_PLOT).d)    )
      ax.set_ylim( 1e-4, 10.0*np.max(Ini_DensProf_dens.in_units(UNIT_D_PLOT).d) )
      ax.set_xscale('log')
      ax.set_yscale('log')

      # set the labels
      ax.legend()
      ax.set_xlabel( r'$r$'+' (%s)'%UNIT_L_PLOT    )
      ax.set_ylabel( r'$\rho$'+' (%s)'%UNIT_D_PLOT )
      ax.set_title( '$t$ = %7.6e code_time'%ds.current_time.in_units('code_time') )

      ax.text( 4.0*np.min(Ini_DensProf_r.in_units(UNIT_L_PLOT).d), 3.2e-3, 'Center Relative Error = %5.1f %%'%(100*(np.abs(prof_dens[field].in_units(UNIT_D_PLOT).d[0]-Ini_DensProf_dens.in_units(UNIT_D_PLOT).d[0])/(Ini_DensProf_dens.in_units(UNIT_D_PLOT).d[0]))), va='center', ha='left' )
      ax.text( 4.0*np.min(Ini_DensProf_r.in_units(UNIT_L_PLOT).d), 8.0e-4, 'Ave( Relative Error ) = %5.1f %%'%(100*np.average(np.abs(prof_dens[field].in_units(UNIT_D_PLOT).d-Ini_DensProf_dens.in_units(UNIT_D_PLOT).d)/(Ini_DensProf_dens.in_units(UNIT_D_PLOT).d+1e-10))), va='center', ha='left' )
      ax.text( 4.0*np.min(Ini_DensProf_r.in_units(UNIT_L_PLOT).d), 2.0e-4, 'Max( Relative Error ) = %5.1f %%'%(100*np.max(    np.abs(prof_dens[field].in_units(UNIT_D_PLOT).d-Ini_DensProf_dens.in_units(UNIT_D_PLOT).d)/(Ini_DensProf_dens.in_units(UNIT_D_PLOT).d+1e-10))), va='center', ha='left' )

      ErrorFile = open( 'DensityProfileErrorEvolution', 'a' )
      ErrorFile.write( '   %9.8e   %9.8e   %9.8e   %9.8e\n'%( ds.current_time.in_units('code_time').d,
                                                              (np.abs(prof_dens[field].in_units(UNIT_D_PLOT).d[0]-Ini_DensProf_dens.in_units(UNIT_D_PLOT).d[0])/(Ini_DensProf_dens.in_units(UNIT_D_PLOT).d[0])),
                                                              np.average(np.abs(prof_dens[field].in_units(UNIT_D_PLOT).d-Ini_DensProf_dens.in_units(UNIT_D_PLOT).d)/(Ini_DensProf_dens.in_units(UNIT_D_PLOT).d+1e-10)),
                                                              np.max(    np.abs(prof_dens[field].in_units(UNIT_D_PLOT).d-Ini_DensProf_dens.in_units(UNIT_D_PLOT).d)/(Ini_DensProf_dens.in_units(UNIT_D_PLOT).d+1e-10)) ) )
      ErrorFile.close()

      # set the grid and ticks
      ax.grid()
      ax.xaxis.set_ticks_position('both')
      ax.yaxis.set_ticks_position('both')
      ax.tick_params( which='both',direction='in' )

      # save the figure
      plt.tight_layout( pad=0.1, h_pad=0.1, w_pad=0.1 )
      fig.savefig( 'fig_%s_Profile_%s.png'%(ds, field), dpi=dpi )
      fig.clear()
