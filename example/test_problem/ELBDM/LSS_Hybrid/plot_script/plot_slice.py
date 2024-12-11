import argparse
import sys
import yt
import numpy as np

# Global parameter settings
fields        = [("gas", "density"), ("gamer", "Phase")]  # Fields to plot
axes          = ["x", "y", "z"]  # Directions
zoom_levels   = [1, 3, 10]  # Zoom levels
max_amr_level = 10  # Maximum level for sampling AMR grid
dpi           = 300
colormap      = 'viridis'


# Load the command-line parameters
parser = argparse.ArgumentParser(description='Plot slices around halo')

# Define command-line arguments
parser.add_argument('-i', action='store', required=False, type=str, dest='prefix',
                    help='path prefix [%(default)s]', default='./')
parser.add_argument('-s', action='store', required=True, type=int, dest='idx_start',
                    help='first data index')
parser.add_argument('-e', action='store', required=True, type=int, dest='idx_end',
                    help='last data index')
parser.add_argument('-d', action='store', required=False, type=int, dest='didx',
                    help='delta data index [%(default)d]', default=1)
parser.add_argument('-m', action='store', required=False, type=int, dest='plot_phase_mod_two_pi',
                    help='plot phase field modulo 2 pi [%(default)d]', default=0)

args = parser.parse_args()

# Print the command-line arguments for reference
print('\nCommand-line arguments:')
print('-------------------------------------------------------------------')
for t in range(len(sys.argv)):
    print(str(sys.argv[t]))
print('')
print('-------------------------------------------------------------------\n')

idx_start             = args.idx_start
idx_end               = args.idx_end
didx                  = args.didx
prefix                = args.prefix
plot_phase_mod_two_pi = args.plot_phase_mod_two_pi

if plot_phase_mod_two_pi:

   def calculate_phase_mod_two_pi(field, data):
      return np.arctan2(np.sin(data["gamer", "Phase"]), np.cos(data["gamer", "Phase"]))

   yt.add_field(("gamer", "PhaseModTwoPi"), function=calculate_phase_mod_two_pi, sampling_type="local", units="")
   fields[1] = ("gamer", "PhaseModTwoPi")

yt.enable_parallelism()

ts = yt.DatasetSeries([prefix + '/Data_%06d' % idx for idx in range(idx_start, idx_end + 1, didx)])

for ds in ts.piter():
    num = '%s' % ds
    num = int(num[5:11])

    ad = ds.all_data()
    ad.max_level = max_amr_level

    loc = ad.quantities.max_location('density')[1:]

    for zoom in zoom_levels:
        for ax in axes:
            for field in fields:
                sz = yt.SlicePlot(ds, ax, field, center=loc, data_source=ad)
                sz.set_cmap(field, colormap)
                if field[1] == "Phase":
                    sz.set_log(("gamer", "Phase"), False)
                if field[1] == "PhaseModTwoPi":
                    sz.set_log(("gamer", "PhaseModTwoPi"), False)
                sz.zoom(zoom)
                sz.set_axes_unit('kpc')
                sz.annotate_timestamp(time_unit='Gyr', redshift=True, corner='upper_right')
                sz.save('Data_%06d_Lv_%02d_Slice_%s_%s_x%d.png' % (num, max_amr_level, ax, field[1], zoom),
                        mpl_kwargs={"dpi": dpi})
                sz.annotate_grids()
                sz.save('Data_%06d_Lv_%02d_Slice_%s_%s_grid_x%d.png' % (num, max_amr_level, ax, field[1], zoom),
                        mpl_kwargs={"dpi": dpi})
