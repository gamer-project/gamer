import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
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
prefix      = args.prefix

colormap    = 'arbre'
dpi         = 150
field       = 'ParDens'


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
for ds in ts.piter():
   p = yt.SlicePlot( ds, 'z', 'ParDens', center='c' )
   field = 'ParDens'
   p.set_background_color( field )
   p.set_unit('ParDens', 'code_mass/(code_length**3)')
   p.set_axes_unit( 'code_length' )
   p.annotate_timestamp( time_unit='code_time', corner='upper_right', text_args={'color':'k'} )
   p.annotate_grids()
   p.save( mpl_kwargs={"dpi":dpi} )
   p.save()
