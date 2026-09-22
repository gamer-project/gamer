import matplotlib.pyplot as plt
import yt


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


# settings for all sub-plots
FileOut = 'fig__profile_merging-cluster.png'

f, ax = plt.subplots( 2, 2 )
#f.subplots_adjust( hspace=0.0, wspace=0.0 )
f.subplots_adjust( wspace=0.5 )
[ f.axes[t].set_xscale( 'log', nonpositive='clip' ) for t in range(0,4,1) ]
[ f.axes[t].set_yscale( 'log', nonpositive='clip' ) for t in range(0,4,1) ]
[ f.axes[t].set_xlim( 5.0e+0, 2.0e+3 ) for t in range(0,4,1) ]
ax[1][0].set_xlabel( r'$r\ [\mathrm{kpc}]$', fontsize=16 )
ax[1][1].set_xlabel( r'$r\ [\mathrm{kpc}]$', fontsize=16 )
for i in range(0,2):
   for j in range(0,2):
      ax[j][i].xaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[2.0,5.0,8.0]) )
      ax[j][i].yaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[2.0,5.0,8.0]) )


# load data
ds_gamer = yt.load( "../Data_000008" )
#ds_gamer = yt.load( "gamer/Data_000100" )
#ds_flash = yt.load( "flash/fiducial_1to1_b0_hdf5_plt_cnt_0100" )


# create a sphere on the max density location
sp_gamer = ds_gamer.sphere( "max", (2.5e3, "kpc") )
#sp_flash = ds_flash.sphere( "max", (2.5e3, "kpc") )


# create a 1D profile object
prof_gamer = yt.create_profile( sp_gamer, 'radius',
                                fields=[ ('gas','density'), ('gas','temperature'), ('gas', 'particle_density_on_grid') ],
                                units = {'radius':'kpc' },
                                weight_field='cell_mass', n_bins=64 )
#prof_flash = yt.create_profile( sp_flash, 'radius',
#                                fields=[ ('gas','density'), ('gas','temperature'), 'pden' ],
#                                units = {'radius':'kpc' },
#                                weight_field='cell_mass', n_bins=64 )


# create arrays to plot
radius_gamer = prof_gamer.x.value
#radius_flash = prof_flash.x.value

dens_gamer = prof_gamer['gas', 'density'    ].value
#dens_flash = prof_flash['gas', 'density'    ].value
temp_gamer = prof_gamer['gas', 'temperature'].value
#temp_flash = prof_flash['gas', 'temperature'].value
pden_gamer = prof_gamer['gas', 'particle_density_on_grid'].value
#pden_flash = prof_flash[       'pden'].value


# remove bins with no data
mask_gamer = temp_gamer > 0.0
#mask_flash = temp_flash > 0.0

radius_gamer = radius_gamer[mask_gamer]
dens_gamer   = dens_gamer  [mask_gamer]
temp_gamer   = temp_gamer  [mask_gamer]
pden_gamer   = pden_gamer  [mask_gamer]

#radius_flash = radius_flash[mask_flash]
#dens_flash   = dens_flash  [mask_flash]
#temp_flash   = temp_flash  [mask_flash]
#pden_flash   = pden_flash  [mask_flash]


# convert gas density to electron number density
dens_gamer *= ne/(amu*mu)
#dens_flash *= ne/(amu*mu)


# convert temperature to keV
temp_gamer *= kB/keV
#temp_flash *= kB/keV


# calculate entropy
entr_gamer = temp_gamer*dens_gamer**(-2.0/3.0)
#entr_flash = temp_flash*dens_flash**(-2.0/3.0)


# convert DM mass density to msun/kpc^3
pden_gamer /= (msun/kpc**3)
#pden_flash /= (msun/kpc**3)


# plot
# electron density
ax[0][0].plot( radius_gamer, dens_gamer, 'b-', lw=1, label='GAMER' )
#ax[0][0].plot( radius_flash, dens_flash, 'r-', lw=1, label='FLASH' )

ax[0][0].set_ylim( 1.0e-5, 2.0e-2 )
ax[0][0].set_ylabel( r'$n_e\ [\mathrm{cm^{-3}}]$', fontsize=16 )

# gas temperature
ax[1][0].plot( radius_gamer, temp_gamer, 'b-', lw=1, label='GAMER' )
#ax[1][0].plot( radius_flash, temp_flash, 'r-', lw=1, label='FLASH' )

ax[1][0].set_ylim( 5.0e+0, 2.0e+1 )
ax[1][0].set_ylabel( r'$T\ [\mathrm{keV}]$', fontsize=16 )

# entropy
ax[0][1].plot( radius_gamer, entr_gamer, 'b-', lw=1, label='GAMER' )
#ax[0][1].plot( radius_flash, entr_flash, 'r-', lw=1, label='FLASH' )

ax[0][1].set_ylim( 1.0e+2, 1.0e+4 )
ax[0][1].set_ylabel( r'$S \ [\mathrm{keV\ cm^2}]$', fontsize=16 )

# dark matter mass density
ax[1][1].plot( radius_gamer, pden_gamer, 'b-', lw=1, label='GAMER' )
#ax[1][1].plot( radius_flash, pden_flash, 'r-', lw=1, label='FLASH' )

ax[1][1].set_ylim( 5.0e+3, 1.0e+8 )
ax[1][1].set_ylabel( r'$\rho_{\mathrm{DM}}\ [\mathrm{M_{\odot}\ kpc^{-3}}]$', fontsize=16 )

# legend
ax[0][0].legend( loc='lower left', numpoints=1, labelspacing=0.2, fontsize=13, handletextpad=0.5,
                 borderpad=0.4, handlelength=1.0 )


# show/save figure
plt.savefig( FileOut, bbox_inches='tight', pad_inches=0.05 )
plt.show()
