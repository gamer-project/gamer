import argparse
import itertools
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Power Spectrum' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='./' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '--compare', action='store_true', dest='compare_with_linear',
                     help='compare with linear evolution (requires data with index 0 in data path) [%(default)d]', default=False )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
    print(str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start           = args.idx_start
idx_end             = args.idx_end
didx                = args.didx
prefix              = args.prefix
compare_with_linear = args.compare_with_linear
dpi                 = 200

if compare_with_linear:
    ds    = yt.load( prefix+'/Data_%06d'%0)
    z0 = ds.current_redshift
    a0 = 1/(1+z0)

    k0, P0 = np.loadtxt( prefix+'/PowerSpec_%06d'%0, skiprows=3, unpack=True)
    ds.close()

yt.enable_parallelism()
ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

zs  = []
ks  = []
Ps  = []
PLs = []
Es  = []

for ds in ts.piter():
    num = '%s'%ds
    num = int(num[5:11])

    z = ds.current_redshift

    k, P = np.loadtxt( prefix+'/PowerSpec_%06d'%num, skiprows=3, unpack=True)
    a = 1/(1 + z)
    zs.append(z)
    ks.append(k)
    Ps.append(P)


    if compare_with_linear:

        PL = (a/a0)**2 * P0
        PLs.append(PL)
        E     = (P - PL) / (PL)
        Es.append(E)

        fig, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1, dpi = dpi, sharex=True)
        plt.suptitle("Power spectra at z = %d" % z)

        ax1.plot(k, k**3*P, label="Simulation" % z)
        ax1.plot(k, k**3*PL, ls = "dashed", lw = 0.5, markersize=3, marker="o", markeredgecolor="k", label="Linear PT" % z)


        ax2.plot(k, E, marker="o", markeredgecolor="k", label="Simulation vs linear PT")

        ax1.legend()
        ax1.set_yscale("log")
        ax1.set_xscale("log")

        ax1.set_ylabel(r'$P(k) \cdot k^3$')

        ax2.set_xscale("log")
        ax2.set_yscale("symlog", linthreshy=1)
        ax2.set_xlabel(r'$k$ in h/Mpc')
        ax2.set_ylabel(r'$(P - P_{lin})/P_{lin}$')
        ax2.legend()
        plt.subplots_adjust(wspace=0, hspace=0)

    else:
        plt.figure(dpi = dpi)
        plt.title("Power spectrum at z = %d" % z )
        plt.plot(k, k**3*P)
        plt.yscale("log")
        plt.xscale("log")
        plt.ylabel(r'$P(k) \cdot k^3$')
        plt.xlabel(r'$k$ in h/Mpc')

    FileOut = 'Data_%06d_PS' %num+'.png'
    plt.savefig(FileOut)
    plt.close()


if compare_with_linear:
    fig, ax = plt.subplots(1, dpi=dpi)
    plt.title("Power spectra (solid) vs linear PT (dashed)")

    for i in range(len(ks)):
        color = next(ax1._get_lines.prop_cycler)['color']
        ax.plot(ks[i], ks[i]**3*Ps[i], color = color, label="z = %d" % zs[i])
        ax.plot(ks[i], ks[i]**3*PLs[i], color = color, lw = 0.5, ls = "dashed", markersize=3, marker="o", markeredgecolor="k")


    ax.legend()
    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.set_ylabel(r'$P(k) \cdot k^3$')
else:
    plt.figure(dpi = dpi)
    plt.title("Power spectra" )
    for i in range(len(ks)):
        plt.plot(ks[i], ks[i]**3*Ps[i], label="z = %d"%zs[i])
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel(r'$P(k) \cdot k^3$')
    plt.xlabel(r'$k$ in h/Mpc')
    plt.legend()

FileOut = 'Data_%d_%d_%d_PS' % (idx_start, idx_end, didx)+'.png'
plt.savefig(FileOut)
plt.close()
