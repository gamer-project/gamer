#====================================================================================================
# Import
#====================================================================================================
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
import MergerUtilities as mus



#====================================================================================================
# Constant
#====================================================================================================
DPI = 300



#====================================================================================================
# Main
#====================================================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description='Plot the data generated from getClusterData.py' )

    parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                         help='path prefix [%(default)s]', default='../' )

    args = parser.parse_args()

    # take note
    print( '\nCommand-line arguments:' )
    print( '-------------------------------------------------------------------' )
    print( ' '.join(sys.argv) )
    print( '-------------------------------------------------------------------\n' )

    prefix = args.prefix


    # S(r)
    data   = np.loadtxt(prefix+"data_S40.txt", skiprows=1)
    Time   = data[:, 0]
    S40    = data[:, 1]
    S_r500 = data[:, 2]

    plt.plot(Time, S40,    label="S(40)",  linewidth=3.0)
    plt.plot(Time, S_r500, label="S_r500", linewidth=3.0)

    plt.ylim(5e-2, 730)
    plt.xlim(0.0, 10.0)
    plt.xlabel("time (Gyr)")
    plt.ylabel(r"$S(r \ [kpc]) \ [\mathrm{keV\ cm^2}]$")
    plt.title("Entropy Evolution")
    plt.legend()
    plt.savefig(prefix+"/S40_evolution.png", dpi=DPI, bbox_inches='tight')


    # Lx(r)
    Time_R, Power = mus.getRecordData(prefix+"/Record__ClusterCenter", [0, 30])
    data   = np.loadtxt(prefix+"/data_Lx.txt", skiprows=1)
    Time   = data[:, 0]
    Lx10   = data[:, 1]
    Lx20   = data[:, 2]
    Lx50   = data[:, 3]
    Lx100  = data[:, 4]
    Lx200  = data[:, 5]
    Lx500  = data[:, 6]
    Lx1000 = data[:, 7]

    plt.plot(Time,   Lx10,   label="Lx(<10 kpc)",   linewidth=1.0)
    plt.plot(Time,   Lx20,   label="Lx(<20 kpc)",   linewidth=1.0)
    plt.plot(Time,   Lx50,   label="Lx(<50 kpc)",   linewidth=1.0)
    plt.plot(Time,   Lx100,  label="Lx(<100 kpc)",  linewidth=1.0)
    plt.plot(Time,   Lx200,  label="Lx(<200 kpc)",  linewidth=1.0)
    plt.plot(Time,   Lx500,  label="Lx(<500 kpc)",  linewidth=1.0)
    plt.plot(Time,   Lx1000, label="Lx(<1000 kpc)", linewidth=1.0)
    plt.plot(Time_R, Power,  label="jet power",     linewidth=1.0)
    plt.ylim(1e41, 1e51)
    plt.xlim(0.0, 10.0)
    plt.yscale("log")
    plt.xlabel("time (Gyr)")
    plt.ylabel("Lx(r) (erg/s)")
    plt.title("Cooling v.s. Heating")
    plt.legend()
    plt.savefig(prefix+"/coolheat_evolution.png", dpi=DPI, bbox_inches='tight')


    # entropy
    data  = np.loadtxt(prefix+"/data_S.txt", skiprows=1)
    Time  = data[:, 0]
    S10   = data[:, 1]
    S20   = data[:, 2]
    S50   = data[:, 3]
    S100  = data[:, 4]
    S200  = data[:, 5]
    S500  = data[:, 6]
    S1000 = data[:, 7]


    plt.plot(Time, S10,   label="S(10)",   linewidth=1.0)
    plt.plot(Time, S20,   label="S(20)",   linewidth=1.0)
    plt.plot(Time, S50,   label="S(50)",   linewidth=1.0)
    plt.plot(Time, S100,  label="S(100)",  linewidth=1.0)
    plt.plot(Time, S200,  label="S(200)",  linewidth=1.0)
    plt.plot(Time, S500,  label="S(500)",  linewidth=1.0)
    plt.plot(Time, S1000, label="S(1000)", linewidth=1.0)

    plt.ylim(5e-2, 3e3)
    plt.xlim(0.0, 10.0)
    plt.yscale("log")
    plt.xlabel("time (Gyr)")
    plt.ylabel("$S(r \ [kpc]) \ [\mathrm{keV\ cm^2}]$")
    plt.title("Entropy Evolution")
    plt.legend()
    plt.savefig(prefix+"/entr_evolution.png", dpi=DPI, bbox_inches='tight')
