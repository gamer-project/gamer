import argparse
import sys
import yt
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.axes_grid1 import AxesGrid

def make_1d_continuous(f):
    for i in range(len(f) - 1):
        while (f[i] - f[i + 1]) > np.pi:
            f[i + 1 :] += 2 * np.pi
        while (f[i] - f[i + 1]) < -np.pi:
            f[i + 1 :] -= 2 * np.pi
    return f


def make_2d_continuous(f):
    for i in range(f.shape[0]):
        make_1d_continuous(f[i, :])
    for i in range(f.shape[1]):
        make_1d_continuous(f[:, i])
    return f


def make_3d_continuous(f):
    for i in range(f.shape[0]):
        make_2d_continuous(f[i, :, :])
    for i in range(f.shape[1]):
        make_2d_continuous(f[:, i, :])
    for i in range(f.shape[2]):
        make_2d_continuous(f[:, :, i])
    return f


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
parser.add_argument( '--plot-grid', dest ="plot_grid", action="store_true")
parser.set_defaults( plot_grid = False)
parser.add_argument( '--use_phase', dest ="use_phase", action="store_true")
parser.set_defaults( use_phase = False)
parser.add_argument( '--3d', dest ="make_3d_plot", action="store_true")
parser.set_defaults( make_3d_plot = False)
parser.add_argument( '--reim', dest ="plot_reim", action="store_true")
parser.set_defaults( plot_reim = False)
###
args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

plot_grid      = args.plot_grid
use_phase      = args.use_phase
make_3d_plot   = args.make_3d_plot
make_reim_plot = args.plot_reim

colormap    = 'viridis'
dpi         = 150

if use_phase:
     def _phase(field, data):
          return data["gamer", "Phase"]
else:
     def _phase(field, data):
          return np.arctan2(data["gamer", "Imag"], data["gamer", "Real"])

     yt.add_field(
          name=("gamer", "Phase"),
          function=_phase,
          sampling_type="local"
     )

for idx in range(idx_start, idx_end+1, didx):
          ds = yt.load("Data_%06d"%idx)
          try:
            ds.force_periodicity()
          except:
            ds.periodicity = (True, True, True)

          grad_fields = ds.add_gradient_fields(("gas", "density"))
          grad_fields = ds.add_gradient_fields(("gamer", "Real"))
          grad_fields = ds.add_gradient_fields(("gamer", "Imag"))



          if make_3d_plot:
                     axes = ["x", "y", "z"]
          else:
                     axes = ["z"]

          for myax in axes:
                     fig = plt.figure(dpi = dpi, figsize=(24, 12))

                     # See https://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
                     # These choices of keyword arguments produce a four panel plot that includes
                     # four narrow colorbars, one for each plot.  Axes labels are only drawn on the
                     # bottom left hand plot to avoid repeating information and make the plot less
                     # cluttered.
                     grid = AxesGrid(
                          fig,
                          (0.075, 0.075, 0.85, 0.85),
                          nrows_ncols=(3, 2),
                          axes_pad=(0.2, 0.0),
                          label_mode="L",
                          share_all=True,
                          cbar_location="right",
                          cbar_mode="edge",
                          direction="row",
                          cbar_size="3%",
                          cbar_pad="0%",
                     )


                     if make_reim_plot:
                        fields = [
                              ("gamer", "Real"),
                              ("gamer", "Real_gradient_magnitude"),
                              ("gamer", "Imag_gradient_magnitude"),
                         ]
                     else:
                         fields = [
                              ("gas", "density"),
                              ("gas", "density_gradient_magnitude"),
                              ("gamer", "Phase"),
                         ]


                     pz = yt.SlicePlot( ds, myax, fields)
                     if not make_reim_plot:
                         pz.set_log(("gamer", "Phase"), False)

                     pz.annotate_grids( periodic=False )

                     for field in fields:
                        pz.set_cmap( field, colormap )

                     # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
                     # axes.
                     for i, field in enumerate(fields):
                          plot = pz.plots[field]
                          plot.figure = fig
                          plot.axes = grid[2*i].axes
                          plot.cax = grid.cbar_axes[i]

                     pz2 = yt.SlicePlot( ds, myax, fields)
                     if not make_reim_plot:
                         pz2.set_log(("gamer", "Phase"), False)

                     for field in fields:
                        pz2.set_cmap( field, colormap )
                     # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
                     # axes.
                     for i, field in enumerate(fields):
                          plot = pz2.plots[field]
                          plot.figure = fig
                          plot.axes = grid[2*i+1].axes

                     # Finally, redraw the plot on the AxesGrid axes.
                     pz._setup_plots()
                     pz2._setup_plots()
                     plt.savefig("Data_%06d_Slice_%s_density.png" % (idx, myax))
                     plt.close()
