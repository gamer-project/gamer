import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Density profile' )

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

field       = 'density'
center_mode = 'max'
dpi         = 150
nbin        = 32

yt.enable_parallelism()
ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   sp = ds.sphere( center_mode, 0.5*ds.domain_width.to_value().max() )

   prof = yt.ProfilePlot( sp, 'radius', field, weight_field='cell_volume', n_bins=nbin )
   prof.set_unit( 'radius', 'kpc' )
   prof.set_xlim( 5.0e-1, 1.0e2 )
#  prof.set_ylim( field, 1.0e-6, 1.0e0 )

   prof.save( mpl_kwargs={"dpi":dpi} )
