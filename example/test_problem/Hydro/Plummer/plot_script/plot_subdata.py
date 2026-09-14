# Plot a SubData_* sub-cadence output (see OPT__OUTPUT_SUBDIV in Input__Parameter):
# gas density slice (GridData) with the massive-particle positions overlaid (Particle).
# SubData_* files load through the standard yt gamer frontend when the tree is output
# (OPT__OUTPUT_SUBDIV_TREE = 1, the default).
#
# Usage: python plot_subdata.py <first SubDumpID> <last SubDumpID> [delta SubDumpID]

import argparse
import yt

parser = argparse.ArgumentParser( description="Plot SubData_* sub-cadence outputs" )
parser.add_argument( "idx_start", type=int, help="first SubDumpID" )
parser.add_argument( "idx_end",   type=int, help="last SubDumpID" )
parser.add_argument( "didx",      type=int, help="delta SubDumpID [1]", nargs="?", default=1 )
args = parser.parse_args()

ts = yt.DatasetSeries( [ "SubData_%06d"%idx
                         for idx in range(args.idx_start, args.idx_end+1, args.didx) ] )

for ds in ts.piter():
   sz = yt.SlicePlot( ds, "z", ("gamer","Dens"), center="c" )
   sz.set_zlim( ("gamer","Dens"), 1.0e-6, 1.0e0 )
   sz.annotate_particles( width=(1.0,"code_length"), p_size=2.0, col="red" )
   sz.annotate_timestamp( time_unit="code_time", corner="upper_right" )
   sz.save( mpl_kwargs={"dpi":150} )
