import argparse
import sys
import yt
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
parser.add_argument( '--plot-grid', dest ="plot_grid", action="store_true")
parser.set_defaults( plot_grid = False)
parser.add_argument( '--plot-phase', dest ="plot_phase", action="store_true")
parser.set_defaults( plot_phase = False)
parser.add_argument( '--plot-3d', dest ="plot_3d", action="store_true")
parser.set_defaults( plot_3d = False)
parser.add_argument( '--plot-reim', dest ="plot_reim", action="store_true")
parser.set_defaults( plot_reim = False)
parser.add_argument( '--plot-grav', dest ="plot_grav", action="store_true")
parser.set_defaults( plot_grav = False)
parser.add_argument( '--comp-phase', dest ="comp_phase", action="store_true")
parser.set_defaults( plot_grav = False)


###
args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
    print( str(sys.argv[t]) ),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start  = args.idx_start
idx_end    = args.idx_end
didx       = args.didx
prefix     = args.prefix

plot_grid  = args.plot_grid
plot_phase = args.plot_phase
plot_3d    = args.plot_3d
plot_reim  = args.plot_reim
plot_grav  = args.plot_grav
comp_phase = args.comp_phase

colormap    = 'arbre'
center_mode = 'c'
dpi         = 150

yt.enable_parallelism()

if comp_phase:
  def _phase(field, data):
          return np.arctan2(data["gamer", "Imag"], data["gamer", "Real"])

  yt.add_field(
      name=("gamer", "Phase"),
      function=_phase,
      sampling_type="local",
      units="",
  )
plot_reim2 = 0

if plot_reim2:
     def _real(field, data):
          return np.cos(data["gamer", "Phase"])

     def _imag(field, data):
          return np.sin(data["gamer", "Phase"])


     yt.add_field(
          name=("gamer", "Imag"),
          function=_imag,
          sampling_type="local",
          units="",
     )

     yt.add_field(
          name=("gamer", "Real"),
          function=_real,
          sampling_type="local",
          units="",
     )

   
for idx in range(idx_start, idx_end+1, didx):
          ds = yt.load("Data_%06d"%idx)# ds = yt.load('t_%06d'%idx)
          ds.force_periodicity()
          grad_fields  = ds.add_gradient_fields(("gas", "density"))
          phase_fields = ds.add_gradient_fields(("gamer", "Real"))
          if plot_grav:
              pot_fields   = ds.add_gradient_fields(("gas", "gravitational_potential"))
          if plot_reim:
              grad_fields = ds.add_gradient_fields(("gamer", "Real"))
              grad_fields = ds.add_gradient_fields(("gamer", "Imag"))



          if plot_3d:
                     axes = ["x", "y", "z"]
          else:
                     axes = ["z"]

          for myax in axes:
                     fig = plt.figure(dpi = 120, figsize=(36 + 12 * plot_grav + 12 * plot_reim, 12))

                     # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
                     # These choices of keyword arguments produce a four panel plot that includes
                     # four narrow colorbars, one for each plot.  Axes labels are only drawn on the
                     # bottom left hand plot to avoid repeating information and make the plot less
                     # cluttered.
                     grid = AxesGrid(
                          fig,
                          (0.075, 0.075, 0.85, 0.85),
                          nrows_ncols=(4+2*plot_grav+4*plot_reim, 2),
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
                          ("gas", "density_gradient_magnitude"),
                          #("gamer", "Phase"),
                          #("gamer", "Phase_gradient_magnitude"),
                         ]
                     
                     if plot_grav:
                         fields.extend([("gas", "gravitational_potential"),
                              ("gas", "gravitational_potential_gradient_magnitude")])
 
                     if plot_reim:
                         fields.extend([
                              ("gamer", "Real"),
                              ("gamer", "Imag")])

                     pz = yt.SlicePlot( ds, myax, fields)
                     pz.set_log(("gamer", "Phase"), False)
                     if plot_reim:
                         pz.set_log(("gamer", "Real"), False)
                         pz.set_log(("gamer", "Imag"), False)

                     pz.annotate_grids( periodic=False )

                     pz.set_cmap( fields, "viridis" )

                     # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
                     # axes.
                     for i, field in enumerate(fields):
                          plot = pz.plots[field]
                          plot.figure = fig
                          plot.axes = grid[2*i].axes
                          plot.cax = grid.cbar_axes[i]

                     pz2 = yt.SlicePlot( ds, myax, fields)
                     pz2.set_log(("gamer", "Phase"), False)
                     if plot_reim:
                         pz2.set_log(("gamer", "Real"), False)
                         pz2.set_log(("gamer", "Imag"), False)

                     pz2.set_cmap( fields, "viridis" )
                     # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
                     # axes.
                     for i, field in enumerate(fields):
                          plot = pz2.plots[field]
                          plot.figure = fig
                          plot.axes = grid[2*i+1].axes

                     # Finally, redraw the plot on the AxesGrid axes.
                     pz._setup_plots()
                     pz2._setup_plots()
                     plt.savefig(f"Data_{idx:06}_density_{myax}_axis.png")

                     plt.close()
