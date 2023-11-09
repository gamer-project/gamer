import argparse
import sys
import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ref: https://yt-project.org/doc/examining/low_level_inspection.html#examining-grid-data-in-a-fixed-resolution-array


# load the command-line parameters
parser = argparse.ArgumentParser( description='Extract power spectrum from GAMER AMR data.' )

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

fields      = [("gas", "density")]   # fields for which to compute power spectrum

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
    num = '%s'%ds
    num = int(num[5:11])
    N   = ds.domain_dimensions * 2 ** lv

    data = ds.smoothed_covering_grid(
        level=lv, left_edge=[0, 0.0, 0.0], dims=N
    )

    for field in fields:
        N = data[field].shape[0]

        d = data[field]

        # compute fft and take average
        power_spectrum = np.mean(np.mean(np.abs(np.fft.fftn(d))**2, axis=1), axis=1)

        # Create a log-log plot of the power spectrum
        plt.figure(figsize=(8, 6))
        k_values = np.fft.fftfreq(N)  # Assuming isotropic grid
        k_values = k_values[:N//2]   # Only use positive k-values
        plt.loglog(k_values, power_spectrum[:N//2], label='Power Spectrum')
        plt.xlabel('k (Wavenumber)')
        plt.ylabel('Power Spectrum')
        plt.title('%s Power Spectrum' % field[1].capitalize())
        plt.grid()
        plt.legend()
        f = (np.char.capitalize(field[1]), num, lv, N)
        plt.savefig('power_spectrum_%s_%06d_lv_%02d_%d.png' % f)
        plt.close()

        # store spectra in array
        np.savetxt('power_spectrum_%s_%06d_lv_%02d_%d.txt' % f, np.column_stack([k_values, power_spectrum[:N//2]]))
