import argparse
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
import numpy as np

# Load the command-line parameters
parser = argparse.ArgumentParser(description='Power Spectrum Analysis')

# Define command-line arguments
parser.add_argument('-i', action='store', required=False, type=str, dest='data_path',
                    help='data path prefix [%(default)s]', default='./')
parser.add_argument('-s', action='store', required=True, type=int, dest='start_index',
                    help='start data index')
parser.add_argument('-e', action='store', required=True, type=int, dest='end_index',
                    help='end data index')
parser.add_argument('-d', action='store', required=False, type=int, dest='index_step',
                    help='index step [%(default)d]', default=1)
parser.add_argument('--compare', action='store_true', dest='compare_with_linear',
                    help='compare with linear evolution (requires data with index 0 in data path) [%(default)d]', default=False)

args = parser.parse_args()

# Print the command-line arguments for reference
print('\nCommand-line arguments:')
print('-------------------------------------------------------------------')
for t in range(len(sys.argv)):
    print(str(sys.argv[t]))
print('')
print('-------------------------------------------------------------------\n')

start_index = args.start_index
end_index = args.end_index
index_step = args.index_step
data_path = args.data_path
compare_with_linear = args.compare_with_linear
dpi = 200

if compare_with_linear:
    ds0 = yt.load(data_path + '/Data_%06d' % 0)
    z0 = ds0.current_redshift
    a0 = 1 / (1 + z0)

    k0, P0 = np.loadtxt(data_path + '/PowerSpec_%06d' % 0, skiprows=3, unpack=True)
    ds0.close()

yt.enable_parallelism()
data_series = yt.DatasetSeries([data_path + '/Data_%06d' % idx for idx in range(start_index, end_index + 1, index_step)])

redshifts = []
wave_numbers = []
power_spectra = []
linear_power_spectra = []
relative_errors = []

for data_set in data_series.piter():
    current_index = int(str(data_set)[5:11])
    current_redshift = data_set.current_redshift

    k, P = np.loadtxt(data_path + '/PowerSpec_%06d' % current_index, skiprows=3, unpack=True)
    scale_factor = 1 / (1 + current_redshift)

    redshifts.append(current_redshift)
    wave_numbers.append(k)
    power_spectra.append(P)

    if compare_with_linear:

        linear_power = (scale_factor / a0) ** 2 * P0
        linear_power_spectra.append(linear_power)
        relative_error = (P - linear_power) / linear_power
        relative_errors.append(relative_error)

        # Create and save comparison plots
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, dpi=dpi, sharex=True)
        plt.suptitle("Power spectra at z = %.2f" % current_redshift)

        ax1.plot(k, k ** 3 * P, label="Simulation")
        ax1.plot(k, k ** 3 * linear_power, ls="dashed", lw=0.5, markersize=3, marker="o", markeredgecolor="k", label="Linear PT")

        ax2.plot(k, relative_error, marker="o", markeredgecolor="k", label="Simulation vs Linear PT")

        ax1.legend()
        ax1.set_yscale("log")
        ax1.set_xscale("log")
        ax1.set_ylabel(r'$P(k) \cdot k^3')

        ax2.set_xscale("log")
        ax2.set_yscale("symlog", linthreshy=1)
        ax2.set_xlabel(r'$k$ in h/Mpc')
        ax2.set_ylabel(r'$(P - P_{lin})/P_{lin}$')
        ax2.legend()
        plt.subplots_adjust(wspace=0, hspace=0)

    else:
        # Create and save power spectrum plots without comparison
        plt.figure(dpi=dpi)
        plt.title("Power spectrum at z = %.2f" % current_redshift)
        plt.plot(k, k ** 3 * P)
        plt.yscale("log")
        plt.xscale("log")
        plt.ylabel(r'$P(k) \cdot k^3$')
        plt.xlabel(r'$k$ in h/Mpc')

    # Save each power spectrum plot
    output_file = 'Data_%06d_PS' % current_index + '.png'
    plt.savefig(output_file)
    plt.close()

if compare_with_linear:
    # Create and save a summary plot comparing power spectra with linear PT
    fig, ax = plt.subplots(1, dpi=dpi)
    plt.title("Power spectra (solid) vs linear PT (dashed)")

    for i in range(len(wave_numbers)):
        color = next(ax._get_lines.prop_cycler)['color']
        ax.plot(wave_numbers[i], wave_numbers[i] ** 3 * power_spectra[i], color=color, label="z = %.2f" % redshifts[i])
        ax.plot(wave_numbers[i], wave_numbers[i] ** 3 * linear_power_spectra[i], color=color, lw=0.5, ls="dashed",
                markersize=3, marker="o", markeredgecolor="k")

    ax.legend()
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel(r'$P(k) \cdot k^3$')
else:
    # Create and save a summary plot of power spectra without comparison
    plt.figure(dpi=dpi)
    plt.title("Power spectra")
    for i in range(len(wave_numbers)):
        plt.plot(wave_numbers[i], wave_numbers[i] ** 3 * power_spectra[i], label="z = %.2f" % redshifts[i])
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel(r'$P(k) \cdot k^3$')
    plt.xlabel(r'$k$ in h/Mpc')
    plt.legend()

output_file = 'Data_%d_%d_%d_PS' % (start_index, end_index, index_step) + '.png'
plt.savefig(output_file)
plt.close()
