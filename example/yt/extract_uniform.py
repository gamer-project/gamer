import argparse
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np

# ref: https://yt-project.org/doc/examining/low_level_inspection.html#examining-grid-data-in-a-fixed-resolution-array


# load the command-line parameters
parser = argparse.ArgumentParser( description='Extract uniform grids from GAMER AMR data.' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                            help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                            help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                            help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False,  type=str, dest='prefix',
                      help='data path prefix [%(default)s]', default='./' )
parser.add_argument( '-lv', action='store', required=False,  type=int, dest='lv',  help='sampling level [%(default)d]', default=0 )



args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print(str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix
lv          = args.lv

fields      = [("gas", "density")]   # fields to extract

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
    num = '%s'%ds
    num = int(num[5:11])
    N   = ds.domain_dimensions * 2 ** lv

    data = ds.smoothed_covering_grid(
        level=lv, left_edge=[0, 0.0, 0.0], dims=N
    )

    for field in fields:
        data[field].tofile('%s_%06d_Lv_%02d_%d_%d_%d.bin'%(np.char.capitalize(field[1]), num, lv, N[0], N[1], N[2]))

