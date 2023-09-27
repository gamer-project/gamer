import argparse
import sys
import yt
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.axes_grid1 import AxesGrid


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot slices of wave function for the ELBDM test' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                            help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                            help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                            help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False,  type=str, dest='prefix',
                      help='data path prefix [%(default)s]', default='./' )


###
args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
    print( str(sys.argv[t]) ),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

colormap  = 'viridis'
dpi       = 150

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
for ds in ts.piter():
     axes = ["z"]

     for myax in axes:
          fig = plt.figure(dpi = dpi, figsize=(24, 12))

          # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
          # These choices of keyword arguments produce a four panel plot that includes
          # four narrow colorbars, one for each plot.  Axes labels are only drawn on the
          # bottom left hand plot to avoid repeating information and make the plot less
          # cluttered.
          grid = AxesGrid(
               fig,
               (0.075, 0.075, 0.85, 0.85),
               nrows_ncols=(2, 2),
               axes_pad=(0.2, 0.0),
               label_mode="L",
               share_all=True,
               cbar_location="right",
               cbar_mode="edge",
               direction="row",
               cbar_size="3%",
               cbar_pad="0%",
          )


          fields = [
               ("gas", "density"),
               ("gamer", "Phase"),
          ]


          pz = yt.SlicePlot( ds, myax, fields)
          pz.set_log(("gamer", "Phase"), False)

          pz.annotate_grids( periodic=False )

          pz.set_cmap( fields, colormap )

          # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
          # axes.
          for i, field in enumerate(fields):
               plot = pz.plots[field]
               plot.figure = fig
               plot.axes = grid[2*i].axes
               plot.cax = grid.cbar_axes[i]

          pz2 = yt.SlicePlot( ds, myax, fields)
          pz2.set_log(("gamer", "Phase"), False)
          pz2.set_cmap( fields, colormap )
          # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
          # axes.
          for i, field in enumerate(fields):
               plot = pz2.plots[field]
               plot.figure = fig
               plot.axes = grid[2*i+1].axes

          # Finally, redraw the plot on the AxesGrid axes.
          pz._setup_plots()
          pz2._setup_plots()

          DumpID = ds.parameters["DumpID"]
          plt.savefig("Data_%06d_density_%s_axis.png" % (DumpID, myax))

          plt.close()
