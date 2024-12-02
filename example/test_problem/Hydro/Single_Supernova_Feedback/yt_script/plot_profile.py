import argparse
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

#parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix_1',
#                     help='path prefix [%(default)s]', default='../' )
#parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix_2',
#                     help='path prefix [%(default)s]', default='../' )
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
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
#prefix_1    = args.prefix_1
#prefix_2    = args.prefix_2

center_mode = [0.50390625, 0.50390625, 0.50390625]
dpi         = 150


yt.enable_parallelism()


for i in range(idx_start, idx_end+1, didx):
   #load data
   ts = yt.load( "/work1/yglee/Software/gamer_new/bin/SSN_energy/Data_%06d"%i )
   es = yt.load( "/work1/yglee/Software/enzo-dev/run/StarParticle/StarParticleSingleTest/DD%04d/data%04d"%(i, i) )

   #add new field for weighting
   def _density_square( field, data ):
      return data["density"]**2

   es.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="g**2/cm**6" )
   ts.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="g**2/cm**6" )


   #create sphere
   GAMER_sphere = ts.sphere( center_mode, 0.5*ts.domain_width.to_value().max() )
   Enzo_sphere  = es.sphere( center_mode, 0.5*es.domain_width.to_value().max() )

   #create profiles

   #density
   GAMER_profile = yt.create_profile( GAMER_sphere, 'radius', fields='density', weight_field='cell_volume', n_bins=32 )
   Enzo_profile  = yt.create_profile( Enzo_sphere,  'radius', fields='density', weight_field='cell_volume', n_bins=32 )

   #temperature
   #GAMER_profile = yt.create_profile( GAMER_sphere, 'radius', fields='temperature', weight_field='density_square', n_bins=32 )
   #Enzo_profile  = yt.create_profile( Enzo_sphere,  'radius', fields='temperature', weight_field='density_square', n_bins=32 )

   # create the figure
   fig = plt.figure()
   ax  = fig.add_subplot(111)

   # plot the profiles

   #density
   ax.plot( GAMER_profile.x.in_units('cm').d, GAMER_profile['density'].in_units('g/cm**3').d, label='GAMER' )
   ax.plot( Enzo_profile.x.in_units('cm').d, Enzo_profile['density'].in_units('g/cm**3').d, label='Enzo' )

   ax.set_xlim( 1.0e19, 1.0e21 )
   ax.set_ylim( 1.0e-27, 1.0e-23 )
   ax.set_xscale( 'log' )
   ax.set_yscale( 'log' )
   ax.legend()
   ax.set_xlabel( 'radius (cm)' )
   ax.set_ylabel( 'density (g/cm**3)' )

   #temperature
   #ax.plot( GAMER_profile.x.in_units('cm').d, GAMER_profile['temperature'].in_units('K').d, label='GAMER' )
   #ax.plot( Enzo_profile.x.in_units('cm').d, Enzo_profile['temperature'].in_units('K').d, label='Enzo' )

   #ax.set_xlim( 1.0e19, 1.0e21 )
   #ax.set_ylim( 1.0e1, 1.0e7 )
   #ax.set_xscale( 'log' )
   #ax.set_yscale( 'log' )
   #ax.legend()
   #ax.set_xlabel( 'radius (cm)' )
   #ax.set_ylabel( 'temperature (K)' )

   # save the figure
   fig.savefig( 'density_profiles_%06d.png'%i )

