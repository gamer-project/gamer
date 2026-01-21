#====================================================================================================
# Import
#====================================================================================================
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import yt
import h5py

sys.dont_write_bytecode = True
import MergerUtilities as mus



#====================================================================================================
# Constant
#====================================================================================================
# code units in CGS
UNIT_L = 3.08567758149000e+24
UNIT_M = 1.98850000000000e+47
UNIT_T = 3.15569252000000e+17
UNIT_V = 9.77813130377481e+06
UNIT_D = 6.76821209430561e-27
UNIT_E = 1.90124167292092e+61
UNIT_P = 6.47121291670164e-13

# some constants in CGS
amu  = 1.660539040e-24
mu   = 0.59242761692650336
ne   = 0.52                # electron number fraction
keV  = 1.6021766208e-9
kB   = 1.38064852e-16
msun = 1.9885e33
kpc  = 3.08567758149e21

G       = 6.6738e-8
yr      = 31556926.0
mp      = 1.672621898e-24
sigma_T = 6.65e-25
c       = 3e10

R200 = 1812.6  # kpc



#====================================================================================================
# Main
#====================================================================================================
if __name__ == "__main__":
    # 1. load the command-line parameters
    parser = argparse.ArgumentParser( description='Generate major cluster data' )

    parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                         help='path prefix [%(default)s]', default='../' )
    parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                         help='first data index' )
    parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                         help='last data index' )
    parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                         help='delta data index [%(default)d]', default=1 )

    args = parser.parse_args()

    # take note
    print( '\nCommand-line arguments:' )
    print( '-------------------------------------------------------------------' )
    print( ' '.join(sys.argv) )
    print( '-------------------------------------------------------------------\n' )

    idx_start = args.idx_start
    idx_end   = args.idx_end
    didx      = args.didx
    prefix    = args.prefix

    # 2. get cluster 0 center
    cluster_id = 0
    cenX, cenY, cenZ = mus.getClusterCen( prefix, cluster_id, idx_start, idx_end, didx )

    # 3. load data
    def _spectral_weighting(field, data):
        return (data["gas", "density"]*ne/(amu*mu))**2*data["gas", "temperature"]**(-0.75)
    yt.add_field(("gas", "spectral_weighting"),function=_spectral_weighting,sampling_type="local",units='g**2/(K**(3/4)*cm**6)')

    yt.enable_parallelism()
    ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

    # 4. prepare arrays
    S40, S_r500 = [], []
    S10, S20, S50, S100, S200, S500, S1000 = [], [], [], [], [], [], []
    Lx10, Lx20, Lx50, Lx100, Lx200, Lx500, Lx1000 = [], [], [], [], [], [], []
    tcool, deltaT = [], []
    Time = []

    # 5. get profile data
    ii = -1+idx_start
    for ds in ts.piter():
        ii += 1
        sp_gamer = ds.sphere( [cenX[ii], cenY[ii], cenZ[ii]], (1.0e3, "kpc") )
        prof_gamer = yt.create_profile( sp_gamer, 'radius',
                                        fields=[ ('gas','density'), ('gas','temperature'), ('gas', 'pressure'), ('gamer', 'TCool') ],
                                        units = {'radius':'kpc'},
                                        weight_field='cell_mass', n_bins=2000 )
        radius_gamer = prof_gamer.x.value
        dens_gamer   = prof_gamer['gas', 'density'    ].value
        temp_gamer   = prof_gamer['gas', 'temperature'].value
        pres_gamer   = prof_gamer['gas', 'pressure'   ].value  # in_units("g/s**2/cm")
        tcool_gamer  = prof_gamer['gamer', 'TCool'    ].value
        mask_gamer   = temp_gamer > 0.0
        radius_gamer = radius_gamer[mask_gamer]
        dens_gamer   = dens_gamer  [mask_gamer]
        temp_gamer   = temp_gamer  [mask_gamer]
        pres_gamer   = pres_gamer  [mask_gamer]
        tcool_gamer  = tcool_gamer [mask_gamer]
        entr = (temp_gamer*kB/keV)*(dens_gamer*ne/(amu*mu))**(-2.0/3.0)
        emis = 1.5*pres_gamer/(tcool_gamer*UNIT_T)   # in (erg/s/cm^3)

        # integrated X-ray emissivity
        N_emis = len(emis)
        Emis = [0.0] * N_emis
        for i in range(N_emis):
            radius    = radius_gamer[:i+1]*kpc   # in cm
            integrand = emis[:i+1]*4*np.pi*radius**2
            Emis[i]   = np.trapz(integrand, radius)

        Time.append(ds.current_time.to('Gyr').value)
        S40.append(np.interp(40, radius_gamer, entr))
        S_r500.append(np.interp(14.8512, radius_gamer, entr))
        Lx10.append(np.interp(10, radius_gamer, Emis))
        Lx20.append(np.interp(20, radius_gamer, Emis))
        Lx50.append(np.interp(50, radius_gamer, Emis))
        Lx100.append(np.interp(100, radius_gamer, Emis))
        Lx200.append(np.interp(200, radius_gamer, Emis))
        Lx500.append(np.interp(500, radius_gamer, Emis))
        Lx1000.append(np.interp(1000, radius_gamer, Emis))
        S10.append(np.interp(10, radius_gamer, entr))
        S20.append(np.interp(20, radius_gamer, entr))
        S50.append(np.interp(50, radius_gamer, entr))
        S100.append(np.interp(100, radius_gamer, entr))
        S200.append(np.interp(200, radius_gamer, entr))
        S500.append(np.interp(500, radius_gamer, entr))
        S1000.append(np.interp(1000, radius_gamer, entr))
        tcool.append(np.interp(40, radius_gamer, tcool_gamer))

        di_1 = ds.disk( [cenX[ii], cenY[ii], cenZ[ii]], [0.0, 0.0, 1.0], (30.0, "kpc"),     (500.0, "kpc") )
        di_2 = ds.disk( [cenX[ii], cenY[ii], cenZ[ii]], [0.0, 0.0, 1.0], (0.1*R200, "kpc"), (500.0, "kpc") )
        di_3 = ds.disk( [cenX[ii], cenY[ii], cenZ[ii]], [0.0, 0.0, 1.0], (0.2*R200, "kpc"), (500.0, "kpc") )
        cy_in  = di_2 - di_1
        cy_out = di_3 - di_2
        prof_in  = yt.create_profile( cy_in , 'radius', fields= ('gas','temperature'), units = {'radius':'kpc' },
                                      weight_field=("gas", "spectral_weighting"), n_bins=1 )
        prof_out = yt.create_profile( cy_out , 'radius', fields= ('gas','temperature'), units = {'radius':'kpc' },
                                      weight_field=("gas", "spectral_weighting"), n_bins=1 )
        temp_in  = prof_in['gas', 'temperature'].value
        temp_out = prof_out['gas', 'temperature'].value
        T_inner = temp_in[0]
        T_outer = temp_out[0]

        deltaT.append((T_outer-T_inner)/T_outer)

    # Save S(r) to file
    data_to_save = np.column_stack((Time, S40, S_r500))
    np.savetxt( prefix+'data_S40.txt', data_to_save, delimiter='\t', header='# Time\tS40\tS_r500', comments='' )

    # Save Lx(r) to file
    data_to_save = np.column_stack((Time, Lx10, Lx20, Lx50, Lx100, Lx200, Lx500, Lx1000))
    np.savetxt( prefix+'data_Lx.txt', data_to_save, delimiter='\t', header='# Time\tLx20\tLx50\tLx100\tLx200\tLx500\tLx1000', comments='' )

    # Save S(r) to file
    data_to_save = np.column_stack((Time, S10, S20, S50, S100, S200, S500, S1000))
    np.savetxt( prefix+'data_S.txt', data_to_save, delimiter='\t', header='# Time\tS20\tS50\tS100\tS200\tS500\tS1000', comments='' )

    # Save S(r) to file
    data_to_save = np.column_stack((Time, tcool, deltaT))
    np.savetxt( perfix+'data_tcool_deltaT.txt', data_to_save, delimiter='\t', header='# Time\ttcool\tdeltaT', comments='' )
