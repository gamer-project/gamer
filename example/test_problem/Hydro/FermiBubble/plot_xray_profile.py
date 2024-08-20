#====================================================================================================
# Imports
#====================================================================================================
import numpy as np
import matplotlib.pyplot as plt
import h5py



#====================================================================================================
# Constants
#====================================================================================================
fontSize        = 16
lineWidth       = 2.0
NormalizedConst = 26e26
backGround      = 4.5



#====================================================================================================
# Main
#====================================================================================================
# hf = h5py.File( "./Projected_FRB_Data_000035_XRay.h5", 'r' )

group1 = hf.get('Map')
group2 = hf.get('Info')

UpperBound_bb     = np.array(group2.get('Keys')).tolist()[0]
UpperBound_ll     = np.array(group2.get('Keys')).tolist()[1]
numAzimuthalAngle = np.array(group2.get('Keys')).tolist()[2]
b_max             = np.array(group2.get('Keys')).tolist()[3]
b_min             = np.array(group2.get('Keys')).tolist()[4]
l_max             = np.array(group2.get('Keys')).tolist()[5]
l_min             = np.array(group2.get('Keys')).tolist()[6]

print("UpperBound_bb     = %e" % UpperBound_bb    )
print("UpperBound_ll     = %e" % UpperBound_ll    )
print("numAzimuthalAngle = %e" % numAzimuthalAngle)
print("b_max             = %e" % b_max            )
print("b_min             = %e" % b_min            )
print("l_max             = %e" % l_max            )
print("l_min             = %e" % l_min            )

observ_p60 = np.loadtxt( './observed_data/xray_profile/p60.dat' , usecols=(0,1), unpack=True )
observ_p50 = np.loadtxt( './observed_data/xray_profile/p50.dat' , usecols=(0,1), unpack=True )
observ_p40 = np.loadtxt( './observed_data/xray_profile/p40.dat' , usecols=(0,1), unpack=True )
observ_m60 = np.loadtxt( './observed_data/xray_profile/m60.dat' , usecols=(0,1), unpack=True )
observ_m50 = np.loadtxt( './observed_data/xray_profile/m50.dat' , usecols=(0,1), unpack=True )
observ_m40 = np.loadtxt( './observed_data/xray_profile/m40.dat' , usecols=(0,1), unpack=True )

b_p60_idx     = int ( ( +60 - b_min ) / ( b_max - b_min ) * UpperBound_bb )
b_p50_idx     = int ( ( +50 - b_min ) / ( b_max - b_min ) * UpperBound_bb )
b_p40_idx     = int ( ( +40 - b_min ) / ( b_max - b_min ) * UpperBound_bb )
b_m60_idx     = int ( ( -60 - b_min ) / ( b_max - b_min ) * UpperBound_bb )
b_m50_idx     = int ( ( -50 - b_min ) / ( b_max - b_min ) * UpperBound_bb )
b_m40_idx     = int ( ( -40 - b_min ) / ( b_max - b_min ) * UpperBound_bb )

ray = np.linspace( l_min, l_max, num=UpperBound_ll )

for angleIdx in range(numAzimuthalAngle):
    NumRow = 2
    NumCol = 3

    fig, ax = plt.subplots( NumRow, NumCol, sharex="col", sharey="row" )
    fig.subplots_adjust( hspace=0.08, wspace=0.08 )
    fig.set_size_inches( 12.0, 6.0 )

    frb = [None]*2
    frb[0] = np.array( group1.get("ProjectedXray_08_keV") )[angleIdx]
    frb[1] = np.array( group1.get("ProjectedXray_15_keV") )[angleIdx]

    profile_p60 = frb[0][b_p60_idx,:]*NormalizedConst+backGround
    profile_p50 = frb[0][b_p50_idx,:]*NormalizedConst+backGround
    profile_p40 = frb[0][b_p40_idx,:]*NormalizedConst+backGround
    profile_m40 = frb[0][b_m40_idx,:]*NormalizedConst+backGround
    profile_m50 = frb[0][b_m50_idx,:]*NormalizedConst+backGround
    profile_m60 = frb[0][b_m60_idx,:]*NormalizedConst+backGround

    profile_p60 = np.roll( profile_p60, int(profile_p60.size/2) )
    profile_p50 = np.roll( profile_p50, int(profile_p50.size/2) )
    profile_p40 = np.roll( profile_p40, int(profile_p40.size/2) )
    profile_m40 = np.roll( profile_m40, int(profile_m40.size/2) )
    profile_m50 = np.roll( profile_m50, int(profile_m50.size/2) )
    profile_m60 = np.roll( profile_m60, int(profile_m60.size/2) )

    profile_p60 = np.flip( profile_p60 )
    profile_p50 = np.flip( profile_p50 )
    profile_p40 = np.flip( profile_p40 )
    profile_m40 = np.flip( profile_m40 )
    profile_m50 = np.flip( profile_m50 )
    profile_m60 = np.flip( profile_m60 )

    ax[0][0].plot( ray, profile_p60, '-' , linewidth=lineWidth, label=r'$+10^{\circ}$', color='r' )
    ax[0][1].plot( ray, profile_p50, '-' , linewidth=lineWidth, label=r'$+20^{\circ}$', color='r' )
    ax[0][2].plot( ray, profile_p40, '-' , linewidth=lineWidth, label=r'$+30^{\circ}$', color='r' )

    ax[1][0].plot( ray, profile_m60, '-' , linewidth=lineWidth, label=r'$-10^{\circ}$', color='r' )
    ax[1][1].plot( ray, profile_m50, '-' , linewidth=lineWidth, label=r'$-20^{\circ}$', color='r' )
    ax[1][2].plot( ray, profile_m40, '-' , linewidth=lineWidth, label=r'$-30^{\circ}$', color='r' )

    ax[0][0].plot( observ_p60[0], observ_p60[1], '-' , linewidth=lineWidth, label=r'$+40^{\circ}$', color='black' )
    ax[0][1].plot( observ_p50[0], observ_p50[1], '-' , linewidth=lineWidth, label=r'$+50^{\circ}$', color='black' )
    ax[0][2].plot( observ_p40[0], observ_p40[1], '-' , linewidth=lineWidth, label=r'$+60^{\circ}$', color='black' )

    ax[1][0].plot( observ_m60[0], observ_m60[1], '-' , linewidth=lineWidth, label=r'$-40^{\circ}$', color='black' )
    ax[1][1].plot( observ_m50[0], observ_m50[1], '-' , linewidth=lineWidth, label=r'$-50^{\circ}$', color='black' )
    ax[1][2].plot( observ_m40[0], observ_m40[1], '-' , linewidth=lineWidth, label=r'$-60^{\circ}$', color='black' )

    observ_l_max = np.amax( observ_p60[0] )
    observ_l_min = np.amin( observ_p60[0] )

    b = [ +60, +50, +40, -60, -50, -40 ]

    for i in range(NumRow):
        for j in range(NumCol):
            ax[i][j].invert_xaxis()

            ax[i][j].set_xlim( observ_l_min, observ_l_max )

            ax[i][j].tick_params( which='both', direction='in' )
            ax[i][j].tick_params( labelsize=fontSize )
            ax[i][j].tick_params( bottom=True, top=True, left=True, right=True )

            ax[i][j].text( 0.78, 0.95, format(b[i*NumCol+j], '+d')+"$^{\circ}$",
                           horizontalalignment='left', verticalalignment='top',
                           bbox=dict(facecolor='blue', alpha=0.2, boxstyle="round", edgecolor='none'),
                           transform=ax[i][j].transAxes, fontsize=fontSize )

            if i == NumRow-1:
                ax[i][j].set_xlabel( r"Galactic longitude ($^{\circ}$)", fontsize=fontSize )

            if j == 0:
                ax[i][j].set_ylabel(r"Photons s$^{-1}$ deg$^{-2}$", fontsize=fontSize)

            if i == 0: ax[i][j].set_ylim( 3, 21 )
            if i == 1: ax[i][j].set_ylim( 3, 11 )

    plt.savefig( "XRay_profile_0.8keV_%03d.png"%angleIdx, bbox_inches='tight', format="png" )
    plt.show()
    plt.close()
