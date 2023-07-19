import argparse
import sys
import yt
import numpy as np


# global parameter settings

fields      = [("gas", "density"), ("gamer", "Phase")]   # fields to plot
axes        = ["x", "y", "z"]                            # directions
zooms       = [1, 3, 10]                                 # zoom levels
lv          = 10                                         # maximum level for sampling AMR grid
dpi         = 300
colormap    = 'algae'


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot slices around halo' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='./' )
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
   print(str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix


def _phasemod2( field, data ):
   return np.arctan2(np.sin(data["gamer", "Phase"]), np.cos(data["gamer", "Phase"]))

yt.add_field( ("gamer", "PhaseMod2"), function=_phasemod2, sampling_type="local", units="")

yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
   num = '%s'%ds
   num = int(num[5:11])

   ad  = ds.all_data()
   ad.max_level = lv

   loc = ad.quantities.max_location('density')[1:]

   for zoom in zooms:
         for ax in axes:
            for field in fields:
               sz = yt.SlicePlot( ds, ax, field, center=loc, data_source=ad )
               sz.set_cmap( field, colormap )
               if field[1] == "Phase":
                  sz.set_log(("gamer", "Phase"), False)
               if field[1] == "PhaseMod2":
                  sz.set_log(("gamer", "PhaseMod2"), False)
               sz.zoom(zoom)
               sz.set_axes_unit( 'kpc' )
               sz.annotate_timestamp( time_unit='Gyr', redshift=True, corner='upper_right' )
               sz.save('Data_%06d_Lv_%02d_Slice_%s_%s_x%d.png'%(num, lv, ax, field[1], zoom), mpl_kwargs={"dpi":dpi} )
               sz.annotate_grids()
               sz.save('Data_%06d_Lv_%02d_Slice_%s_%s_grid_x%d.png'%(num, lv, ax, field[1], zoom), mpl_kwargs={"dpi":dpi} )