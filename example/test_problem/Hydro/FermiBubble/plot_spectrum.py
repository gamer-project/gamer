#====================================================================================================
# Imports
#====================================================================================================
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import h5py
import matplotlib.font_manager as font_manager
import os.path



#====================================================================================================
# Constants
#====================================================================================================
Hz2GHz             = 1e9
CGS2JANSKY         = 1e-23
DeltaL             = 1.687e21
eV2GeV             = 1e9
CRIndex            = 2.6
CRIndices          = [ 2.2, 2.4, 2.6 ]

TickLabelSize      = 35
LabelSize          = 35
LegendSize         = 24
Pad                = 10
MajorTickLength    = 16
MinorTickLength    = 8

linestyle          = [ "-", "--", "-.", ":" ]
linecolor          = [ "red", "orange", "green", "blue" ]

font               = font_manager.FontProperties( family='monospace', style='normal', size=LegendSize )

path_gamma_ray     = "gamma_ray_spectrum/"
path_synchrotron   = "synchrotron_spectrum/"



#====================================================================================================
# Functions
#====================================================================================================
def trans_b( b, inverse=False ):
    if inverse:
        # transform from [0, 360] to [+90, -90]
        b = ( 180 - b ) * 0.5
    else:
        # transrorm from [+90, -90] to [0, 360]
        b = 180 - 2*b
    return b

def trans_l( l, inverse=False ):
    if inverse:
        # transform from [0, 360] to [-180, +180]
        l = l - 180
    else:
        # transform from [-180, +180] to [0, 360]
        l = 180 + l
    return l

def CRIndex2ScaleCREngy( CRIndex ):
    if   CRIndex == 2.8: ScaleCREngy = 8000
    elif CRIndex == 2.6: ScaleCREngy = 58
    elif CRIndex == 2.4: ScaleCREngy = 7
    elif CRIndex == 2.2: ScaleCREngy = 1.18
    elif CRIndex == 2.0: ScaleCREngy = 1
    else: raise ValueError("Please assign a ScaleCREngy !!")
    return ScaleCREngy


gammaray_B_max = np.array( [50, 30, -10, -30] )
gammaray_B_min = np.array( [30, 10, -30, -50] )

synchrotron_B_max = np.array( [+30, -20] )
synchrotron_B_min = np.array( [+20, -30] )

L_max = +10
L_min = -10

gammaray_B_max = trans_b( gammaray_B_max )
gammaray_B_min = trans_b( gammaray_B_min )

synchrotron_B_max = trans_b( synchrotron_B_max )
synchrotron_B_min = trans_b( synchrotron_B_min )

L_max = trans_l( L_max )
L_min = trans_l( L_min )


fig, ax = plt.subplots( 3, 2 )

fig.set_size_inches( 28, 20 )
fig.subplots_adjust( wspace=0.3, hspace=0.05 )

Emin, Emid, Emax = np.loadtxt( "energyBinsTable", comments='#', usecols=(0,1,2), unpack=True )

###################################
########### Gamma-ray #############
###################################
EngyUpper, spectrumUpper = np.loadtxt( "observedEnergySpectrum-Upper", comments='#', usecols=(0,1), unpack=True, delimiter=',' )
EngyLower, spectrumLower = np.loadtxt( "observedEnergySpectrum-Lower", comments='#', usecols=(0,1), unpack=True, delimiter=',' )

Engy             = np.concatenate( (EngyLower[::-1],     EngyUpper    ) )
ObservedSpectrum = np.concatenate( (spectrumLower[::-1], spectrumUpper) )

for i in range(len(CRIndices)):
    spectrumMultipleB = []
    EmidMultiple      = []

    for b_min, b_max in zip(gammaray_B_max, gammaray_B_min):
        spectrumOneB  = np.array([])
        EmidOneB      = np.array([])

        for binIdx in range(len(Emid)):
            binWidth = Emax[binIdx] - Emin[binIdx]
            filename = path_gamma_ray + "Projected_FRB_Data_000035_%9.4e_1.1e6_"%Emid[binIdx] + str(CRIndices[i]) + ".h5"
            if not os.path.isfile( filename ): continue
            fs = h5py.File( filename, "r" )
            group = fs.get('Map')
            gmap = np.array( group.get('ProjectedLeptonicGammaRay') )
            gmap = gmap[0,:,:]
            gmap = np.roll( gmap, int(gmap.shape[1]/2), axis=1 )
            spectrumOneB = np.append( spectrumOneB, CRIndex2ScaleCREngy(CRIndices[i]) * DeltaL * np.average( gmap[b_min:b_max,L_min:L_max] * binWidth / eV2GeV ) )
            EmidOneB = np.append( EmidOneB, Emid[binIdx] )

        spectrumMultipleB.append( spectrumOneB )
        EmidMultiple.append( EmidOneB )

    for b_idx in range(len(gammaray_B_max)):
       label = r"$%+d^{\circ}<b<%+d^{\circ}$, |l|<10$^{\circ}$"%( trans_b(gammaray_B_min[b_idx], True), trans_b(gammaray_B_max[b_idx], True) )
       ax[i][1].plot( EmidMultiple[b_idx]/eV2GeV, spectrumMultipleB[b_idx], label=label, lw=5, ls=linestyle[b_idx], color=linecolor[b_idx] )

    ax[i][1].fill( Engy, ObservedSpectrum, "silver" )

for i in range(3):
    ax[i][1].set_ylim( [2e-8, 2.5e-6] )
    ax[i][1].set_xlim( [0.1, 700.0] )
    ax[i][1].set_ylabel( "$E^{2}dN/dE$\n(GeV cm$^{-2}$s$^{-1}$sr$^{-1}$)", fontsize=LabelSize )

ax[2][1].set_xlabel( "Photon energy (GeV)", fontsize=LabelSize )

###################################
########### Synchrotron ###########
###################################
WMAPData_x, WMAPData_y = np.loadtxt( "WMAPData-1", comments='#', usecols=(0,1), unpack=True, delimiter=',' )
Engy, ObservedSpectrum = np.loadtxt( "WMAPData-2", comments='#', usecols=(0,1), unpack=True, delimiter=',' )

for i in range(len(CRIndices)):
    spectrumMultipleB = []
    EmidMultiple      = []

    for b_min, b_max in zip(synchrotron_B_max, synchrotron_B_min):
        spectrumOneB = np.array([])
        EmidOneB     = np.array([])

        for binIdx in range(len(Emid)):
            binWidth = Emax[binIdx] - Emin[binIdx]
            filename = path_synchrotron + "Projected_FRB_Data_000035_Synchrotron_%9.4e_1.1e6_"%Emid[binIdx] + str(CRIndices[i]) + ".h5"

            if not os.path.isfile( filename ): continue

            fs = h5py.File( filename, "r" )
            group = fs.get('Map')
            gmap = np.array( group.get('ProjectedSynchrotron') )
            gmap = gmap[0,:,:]
            gmap = np.roll( gmap, int(gmap.shape[1]/2), axis=1 )
            spectrumOneB = np.append( spectrumOneB, CRIndex2ScaleCREngy(CRIndices[i])*DeltaL/1e40*np.average( gmap[b_min:b_max,L_min:L_max] ) / (1e3*CGS2JANSKY) )
            EmidOneB = np.append( EmidOneB, Emid[binIdx] )

        spectrumMultipleB.append( spectrumOneB )
        EmidMultiple.append( EmidOneB )

    for b_idx in range(len(synchrotron_B_max)):
        label = r"$%+d^{\circ}<b<%+d^{\circ}$, |l|<10$^{\circ}$"%( trans_b(synchrotron_B_min[b_idx], True), trans_b(synchrotron_B_max[b_idx], True) )
        ax[i][0].plot( EmidMultiple[b_idx]/Hz2GHz, spectrumMultipleB[b_idx], label=label, lw=5, ls=linestyle[b_idx], color=linecolor[b_idx] )

for i in range(3):
    ax[i][0].set_ylim( [0.08, 13] )
    ax[i][0].set_xlim( [1, max(EmidMultiple[0]/Hz2GHz)] )
    ax[i][0].set_ylabel( r"$I_{\nu}$ (kJy sr$^{-1}$)", fontsize=LabelSize )
    ax[i][0].fill( Engy, ObservedSpectrum, "silver" )
    yerr = np.abs( WMAPData_y[0] - np.average(WMAPData_y) )
    ax[i][0].errorbar( WMAPData_x[0], np.average(WMAPData_y), yerr=yerr, label='WMAP', capsize=10, elinewidth=5, capthick=5 )

ax[2][0].set_xlabel( "Photon frequency (GHz)", fontsize=LabelSize )

for i in range(3):
    for j in range(2):
        ax[i][j].tick_params( axis='both', which='major', colors='k', direction='in', top=True, bottom=True, left=True, right=True, labelsize=TickLabelSize, pad=Pad, length=MajorTickLength )
        ax[i][j].tick_params( axis='both', which='minor', colors='k', direction='in', top=True, bottom=True, left=True, right=True, labelsize=TickLabelSize, pad=Pad, length=MinorTickLength )
        ax[i][j].set_xscale( 'log' )
        ax[i][j].set_yscale( 'log' )

ax[2][0].legend( loc="lower left",  prop=font, handlelength=2, borderaxespad=1.0 )
ax[2][1].legend( loc="center left", prop=font, handlelength=2, borderaxespad=1.0, bbox_to_anchor=(0.06,0.29) )

ax[0][0].get_xaxis().set_ticks([])
ax[1][0].get_xaxis().set_ticks([])
ax[0][1].get_xaxis().set_ticks([])
ax[1][1].get_xaxis().set_ticks([])

plt.savefig( "Spectrum.png", bbox_inches='tight' )
