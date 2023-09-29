from __future__ import print_function, division, absolute_import

import argparse
import sys
import yt
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting

import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.axes_grid1 import AxesGrid

# Load the command-line parameters
parser = argparse.ArgumentParser(description='Plot slices of wave function for the ELBDM test')

# Define command-line arguments
parser.add_argument('-s', action='store', required=True, type=int, dest='idx_start',
                    help='first data index')
parser.add_argument('-e', action='store', required=True, type=int, dest='idx_end',
                    help='last data index')
parser.add_argument('-d', action='store', required=False, type=int, dest='didx',
                    help='delta data index [%(default)d]', default=1)
parser.add_argument('-i', action='store', required=False, type=str, dest='prefix',
                    help='data path prefix [%(default)s]', default='./')
parser.add_argument('-c', action='store', required=False, type=int, dest='convert_reim_to_phase',
                    help='convert real and imaginary part to phase for wave scheme [%(default)d]', default=0)

args = parser.parse_args()  # Parse the command-line arguments

# Print the command-line arguments for reference
print('\nCommand-line arguments:')
print('-------------------------------------------------------------------')
for t in range(len(sys.argv)):
    print(str(sys.argv[t]), end=' ')
print('')
print('-------------------------------------------------------------------\n')

# Extract arguments from the parsed command-line arguments
idx_start             = args.idx_start
idx_end               = args.idx_end
didx                  = args.didx
prefix                = args.prefix
convert_reim_to_phase = args.convert_reim_to_phase

if convert_reim_to_phase:

   def reim2phase(field, data):
      return np.arctan2(data["gamer", "Imag"], data["gamer", "Real"])

   yt.add_field(("gamer", "Phase"), function=reim2phase, sampling_type="local", units="")

colormap = 'viridis'  # Define the colormap for the plots
dpi = 150  # Define the DPI (dots per inch) for the saved plots

# Create a series of datasets based on data files with indices in the specified range
dataset_series = yt.DatasetSeries([prefix + '/Data_%06d' % idx for idx in range(idx_start, idx_end + 1, didx)])

# Loop through each dataset in the series
for dataset in dataset_series.piter():
     axes = ["z"]  # Specify the axes for slicing (e.g., "z" for z-axis slices)

     for current_axis in axes:
          # Create a new figure for the current slice
          fig = plt.figure(dpi=dpi, figsize=(24, 12))

          # Create a grid of axes for multiple plots
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

          # Define the fields to plot
          fields_to_plot = [
               ("gas", "density"),
               ("gamer", "Phase"),
          ]

          # Create a slice plot for the current dataset and field
          slice_plot = yt.SlicePlot(dataset, current_axis, fields_to_plot)
          slice_plot.set_log(("gamer", "Phase"), False)

          slice_plot.annotate_grids(periodic=False)

          for field in fields_to_plot:
               slice_plot.set_cmap(field, colormap)

          # For each plotted field, associate it with the corresponding AxesGrid axes
          for i, field in enumerate(fields_to_plot):
               plot = slice_plot.plots[field]
               plot.figure = fig
               plot.axes = grid[2 * i].axes
               plot.cax = grid.cbar_axes[i]

          # Create a second slice plot for comparison
          slice_plot_2 = yt.SlicePlot(dataset, current_axis, fields_to_plot)
          slice_plot_2.set_log(("gamer", "Phase"), False)

          for field in fields_to_plot:
               slice_plot_2.set_cmap(field, colormap)

          # Associate the second slice plot with the AxesGrid axes
          for i, field in enumerate(fields_to_plot):
               plot = slice_plot_2.plots[field]
               plot.figure = fig
               plot.axes = grid[2 * i + 1].axes

          # Redraw the plot on the AxesGrid axes
          slice_plot._setup_plots()
          slice_plot_2._setup_plots()

          # Get the DumpID from dataset parameters and save the plot
          dump_id = dataset.parameters["DumpID"]
          plt.savefig("Data_%06d_%s_axis.png" % (dump_id, current_axis))

          # Close the current figure to release resources
          plt.close()
