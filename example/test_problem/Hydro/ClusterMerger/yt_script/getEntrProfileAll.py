#====================================================================================================
# Import
#====================================================================================================
import yt
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
import MergerUtilities as mus



#====================================================================================================
# Constant
#====================================================================================================
mu  = 0.59242761692650336
ne  = 0.52                      # electron number fraction
kB  = yt.units.kb
amu = yt.units.amu
Mpc = yt.units.Mpc



#====================================================================================================
# Function
#====================================================================================================
def _entropy(field, data):
    return kB*data["gas", "temperature"]*(data["gas", "density"]*ne/(amu*mu))**(-2.0/3.0)

def save_prof_to_file(radius, entropy, filename):
    with open(filename, 'w') as file:
        file.write('# radius (kpc) initial_entr (keV*cm^2) m_entr (keV*cm^2) mc_entr (keV*cm^2) mcf_entr (keV*cm^2)\n')
        for r, e_init, e_m, e_mc, e_mcf in zip(radius, entropy['initial'], entropy['m'], entropy['mc'], entropy['mcf']):
            file.write(f"{r} {e_init} {e_m} {e_mc} {e_mcf}\n")

def getEntropyProfile(data, center, target_radius):
    ds   = yt.load(data)
    sp   = ds.sphere(center, (1.0e4, "kpc"))
    prof = yt.create_profile( sp, 'radius', fields=[('gas','entropy'), ('gas','temperature')],
                              units={'radius':'kpc'}, weight_field='cell_mass', n_bins=100 )
    radius = prof.x.value
    entr   = prof['gas','entropy'].value
    temp   = prof['gas','temperature'].value
    mask   = temp > 0.0

    return np.interp(target_radius, radius[mask], entr[mask])



#====================================================================================================
# Main
#====================================================================================================
if __name__ == "__main__":
    yt.add_field( ("gas", "entropy"), function=_entropy, sampling_type="local", units="keV*cm**2" )

    parser = argparse.ArgumentParser( description='Get the cluster entropy profile from `m`, mc`, and `mcf` cases ' )

    parser.add_argument( '-p1', action='store', required=True, type=str, dest='prefix_m',
                         help='`m` simulation path prefix [%(default)s]' )
    parser.add_argument( '-p2', action='store', required=True, type=str, dest='prefix_mc',
                         help='`mc` simulation path prefix [%(default)s]' )
    parser.add_argument( '-p3', action='store', required=True, type=str, dest='prefix_mcf',
                         help='`mcf` simulation path prefix [%(default)s]' )

    args = parser.parse_args()

    # take note
    print( '\nCommand-line arguments:' )
    print( '-------------------------------------------------------------------' )
    print( ' '.join(sys.argv) )
    print( '-------------------------------------------------------------------\n' )

    prefix_m   = args.prefix_m
    prefix_mc  = args.prefix_mc
    prefix_mcf = args.prefix_mcf


    labels = [ 'initial', 'm', 'mc', 'mcf' ]

    # Define a common radius range for interpolation
    common_radius = np.logspace( np.log10(1e0), np.log10(1e4), 100 )
    entropy = {label:np.zeros_like(common_radius) for label in labels}

    # Load initial profile
    cenX, cenY, cenZ = mus.getClusterCen(prefix_m, 0, 0, 100, 1)
    center = [ cenX[0], cenY[0], cenZ[0] ]
    entropy['initial'] = getEntropyProfile( prefix_m+"/Data_000000", center, common_radius)

    # Load profiles for m, mc, mcf
    last_idx = 100
    cenX, cenY, cenZ = mus.getClusterCen(prefix_m, 0, 0, 100, 1)
    center = [ cenX[last_idx], cenY[last_idx], cenZ[last_idx] ]
    entropy['m'] = getEntropyProfile( prefix_m+"/Data_%06d"%last_idx, center, common_radius)

    cenX, cenY, cenZ = mus.getClusterCen(prefix_mc, 0, 0, 100, 1)
    center = [ cenX[last_idx], cenY[last_idx], cenZ[last_idx] ]
    entropy['mc'] = getEntropyProfile( prefix_mc+"/Data_%06d"%last_idx, center, common_radius)

    cenX, cenY, cenZ = mus.getClusterCen(prefix_mcf, 0, 0, 100, 1)
    center = [ cenX[last_idx], cenY[last_idx], cenZ[last_idx] ]
    entropy['mcf'] = getEntropyProfile( prefix_mcf+"/Data_%06d"%last_idx, center, common_radius)

    save_prof_to_file(common_radius, entropy, prefix_mcf+'/entr_prof.txt')
