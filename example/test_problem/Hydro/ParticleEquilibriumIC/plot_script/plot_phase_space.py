import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='PhasePlot of particles' )

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


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

dpi       = 150

yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
for ds in ts.piter():

    for direction in ['y', 'z']:
        p = yt.ParticlePhasePlot( ds, "particle_velocity_x", "particle_velocity_"+direction, ["particle_mass"] )

        # setting for the figure
        p.set_unit( 'particle_velocity_x',          'code_length/code_time' )
        p.set_unit( 'particle_velocity_'+direction, 'code_length/code_time' )
        p.set_unit( 'particle_mass',                'code_mass'             )
        p.set_log(  'particle_mass',                True                    )
        p.set_xlim( -1.0, 1.0 )
        p.set_ylim( -1.0, 1.0 )
        p.set_zlim( 'particle_mass', 1.0e-7, 1.0e-5 )
        p.set_cmap( 'particle_mass', 'viridis'      )

        # save the figure
        p.save( 'fig_%s_PhasePlotV%sVx.png'%(ds,direction), mpl_kwargs={'dpi':dpi} )

    for direction in ['x', 'y', 'z']:
        p = yt.ParticlePhasePlot( ds, "particle_position_x", "particle_velocity_"+direction, ["particle_mass"] )

        # setting for the figure
        p.set_unit( 'particle_position_x',          'code_length'           )
        p.set_unit( 'particle_velocity_'+direction, 'code_length/code_time' )
        p.set_unit( 'particle_mass',                'code_mass'             )
        p.set_log(  'particle_mass',                True                    )
        p.set_xlim(  0.0, 3.0 )
        p.set_ylim( -1.0, 1.0 )
        p.set_zlim( 'particle_mass', 1.0e-7, 1.0e-5 )
        p.set_cmap( 'particle_mass', 'viridis'      )

        # save the figure
        p.save( 'fig_%s_PhasePlotV%sX.png'%(ds,direction), mpl_kwargs={'dpi':dpi} )
