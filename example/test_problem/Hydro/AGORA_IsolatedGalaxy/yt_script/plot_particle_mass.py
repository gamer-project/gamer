import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the particle projection' )

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
prefix      = '../'

colormap    = 'algae'
width_kpc   = 30
depth_kpc   = 30
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
#ts = yt.load( 'Data_??????' )


# define the particle filter for the newly formed stars
def new_star( pfilter, data ):
   filter = data[ "all", "ParCreTime" ] > 0
   return filter

yt.add_particle_filter( "new_star", function=new_star, filtered_type="all", requires=["ParCreTime"] )

AllPar = ( 'all',      'particle_mass' )
NewPar = ( 'new_star', 'particle_mass' )


for ds in ts.piter():

#  add the particle filter
   ds.add_particle_filter( "new_star" )


#  face-on (all particles)
   pz = yt.ParticleProjectionPlot( ds, 'z', AllPar, center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   pz.set_unit( AllPar, 'Msun' )
   pz.set_zlim( AllPar, 1.0e6, 1.0e8 )
   pz.set_cmap( AllPar, colormap )
   pz.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   pz.save( name=ds.basename+'_AllPar', mpl_kwargs={"dpi":dpi} )


#  face-on (new particles)
   pz = yt.ParticleProjectionPlot( ds, 'z', NewPar, center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   pz.set_unit( NewPar, 'Msun' )
   pz.set_zlim( NewPar, 5.0e4, 1.0e6 )
   pz.set_cmap( NewPar, colormap )
   pz.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   pz.save( name=ds.basename+'_NewPar', mpl_kwargs={"dpi":dpi} )


#  edge-on (all particles)
   px = yt.ParticleProjectionPlot( ds, 'x', AllPar, center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   px.set_unit( AllPar, 'Msun' )
   px.set_zlim( AllPar, 1.0e6, 1.0e8 )
   px.set_cmap( AllPar, colormap )
   px.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   px.save( name=ds.basename+'_AllPar', mpl_kwargs={"dpi":dpi} )


#  edge-on (new particles)
   px = yt.ParticleProjectionPlot( ds, 'x', NewPar, center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   px.set_unit( NewPar, 'Msun' )
   px.set_zlim( NewPar, 5.0e4, 1.0e6 )
   px.set_cmap( NewPar, colormap )
   px.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   px.save( name=ds.basename+'_NewPar', mpl_kwargs={"dpi":dpi} )
