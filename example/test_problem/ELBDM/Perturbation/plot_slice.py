import argparse
import sys
import yt

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

args = parser.parse_args()

# Print the command-line arguments for reference
print('\nCommand-line arguments:')
print('-------------------------------------------------------------------')
for t in range(len(sys.argv)):
    print(str(sys.argv[t]), end=' ')
print('')
print('-------------------------------------------------------------------\n')

# Extract command-line arguments
idx_start = args.idx_start
idx_end = args.idx_end
didx = args.didx
prefix = args.prefix

# Constants
colormap = 'arbre'
field = 'Dens'
center_mode = 'c'
dpi = 150

# Enable parallelism
yt.enable_parallelism()

# Load dataset series
ts = yt.load([prefix + '/Data_%06d' % idx for idx in range(idx_start, idx_end + 1, didx)])

for ds in ts.piter():
   # Create a slice plot
   sz = yt.SlicePlot(ds, 'z', field, center=center_mode)
   # Configure plot settings
   sz.set_log(field, False)
   sz.set_cmap(field, colormap)
   sz.annotate_timestamp(corner='upper_right', time_format='t = {time:.4f} {units}')
   # Save the plot
   sz.save(mpl_kwargs={"dpi": dpi})
