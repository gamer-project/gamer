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
parser.add_argument( '-L_box', action='store', required=True,  type=float, dest='L_box',  help='box size in [Mpc/h]' )


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
L_box       = args.L_box

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

        # Compute the 3D Fourier transform of the density data
        density_fft = np.fft.rfftn(d)


        # Compute the power spectrum
        power_spectrum = np.abs(density_fft)**2

        # Calculate the k values, with rfftn last axis undergoes real transform
        kx_lin, ky_lin, kz_lin = np.fft.fftfreq(d.shape[0]), np.fft.fftfreq(d.shape[1]), np.fft.rfftfreq(d.shape[2])

        # Create a grid of k values
        kx, ky, kz = np.meshgrid(kx_lin, ky_lin, kz_lin, indexing='ij')
        k          = np.sqrt(kx**2 + ky**2 + kz**2)


        bins = kz_lin
        bins += (bins[1] - bins[0])/2

        counts, bin_edges = np.histogram(k, bins)
        totals, bin_edges = np.histogram(k, bins, weights=power_spectrum)

        hist     = totals/counts
        volume   = L_box ** 3
        dc_mode  = np.abs(hist[0])
        ave_dens = np.sqrt(dc_mode) / volume
        Norm     = volume / (ave_dens**2 * volume**2)
        spectrum = hist * Norm

        ks       = (bins + 0.5 * (bins[1] - bins[0])) * 2 * np.pi / L_box * N
        ks       = ks[:-1]

        # Create a log-log plot of the power spectrum
        plt.figure(figsize=(8, 6))
        plt.loglog(ks, spectrum, label='Power Spectrum')
        plt.xlabel('k in [h/Mpc]')
        plt.ylabel('P(k)')
        plt.title('%s Power Spectrum' % field[1].capitalize())
        plt.grid()
        plt.legend()
        f = (np.char.capitalize(field[1]), num, lv, N)
        plt.savefig('power_spectrum_%s_%06d_lv_%02d_%d.png' % f)
        plt.close()

        # store spectra in array
        to_save = np.hstack([ks, spectrum])
        np.savetxt('power_spectrum_%s_%06d_lv_%02d_%d.txt' % f, to_save)

