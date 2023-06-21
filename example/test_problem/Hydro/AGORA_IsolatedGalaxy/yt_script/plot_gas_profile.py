import yt
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot various gas profiles' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = '../'
fileout     = 'fig__gas_profile'

disk_normal = [0.0, 0.0, 1.0]
width_kpc   = 30
nbin        = 50
markersize  = 4.0
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
#ts = yt.load( 'Data_??????' )


# loop over all datasets
for ds in ts.piter():

#  calculate various profiles
#  ==================================================================================
#  define center as the location of peak gas density within 1 kpc from the center of gas mass
   v, cen1 = ds.h.find_max( ("gas", "density") )
   sp1  = ds.sphere( cen1, (30.0, "kpc") )
   cen2 = sp1.quantities.center_of_mass( use_gas=True, use_particles=False ).in_units( "kpc" )
   sp2  = ds.sphere( cen2, (1.0, "kpc") )
   cen3 = sp2.quantities.max_location( ("gas", "density") )
   cen  = ds.arr( [cen3[1].d, cen3[2].d, cen3[3].d], 'code_length' )


#  only include the data within a sphere with a radius of width_kpc
   sp = ds.sphere( cen, (0.5*width_kpc, "kpc") )
   sp.set_field_parameter( "normal", disk_normal )


#  (1) gas surface density
   prof     = yt.ProfilePlot( sp, ("index", "cylindrical_r"), ("gas", "cell_mass"), weight_field=None,
                              n_bins=nbin, x_log=False, accumulation=False )
   gas_dens = prof.profiles[0]["cell_mass"].in_units("Msun").d
   radius   = prof.profiles[0].x.in_units("kpc").d

#  convert mass to surface density in Msun/pc^2
   dr = radius[1] - radius[0] # assuming linear bin
   for b in range( len(radius) ):
      area         = np.pi*( (radius[b]+0.5*dr)**2 - (radius[b]-0.5*dr)**2 )
      gas_dens[b] /= area*1.0e6


#  (2) gas temperature
#  define density^2 for calculating the weighted temperature
   def _density_square( field, data ):
      return data["density"]**2
   ds.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="g**2/cm**6" )

   prof     = yt.ProfilePlot( sp, ("index", "cylindrical_r"), ("gas", "temperature"), weight_field="density_square",
                               n_bins=nbin, x_log=False, accumulation=False )
   gas_temp = prof.profiles[0]["temperature"].in_units("K").d


#  (3) gas rotational velocity
#  consider only dense enough gas in order to exclude the gaseous halo
#  --> follow the AGORA analysis script: https://bitbucket.org/mornkr/agora-analysis-script/
   sp_dense = sp.cut_region( ["obj['gas', 'density'].in_units('g/cm**3') > 1.e-25"] )
   prof     = yt.ProfilePlot( sp_dense, ("index", "cylindrical_r"),  ("gas", "cylindrical_tangential_velocity"),
                              weight_field=("gas", "cell_mass"), n_bins=nbin, x_log=False )
   gas_vrot = prof.profiles[0]["cylindrical_tangential_velocity"].in_units("km/s").d


#  (4) gas velocity dispersion
#  --> follow the AGORA analysis script: https://bitbucket.org/mornkr/agora-analysis-script/
   def _local_rotational_velocity_x( field, data ):
      vx = np.zeros( data[("gas", "velocity_x")].shape )
      for r, vrot in zip(radius, gas_vrot):
         idx = np.where( (data[("index", "cylindrical_r")].in_units("kpc") >= (r - 0.5*dr)) &
                         (data[("index", "cylindrical_r")].in_units("kpc") <  (r + 0.5*dr)) )
         vx[idx] = -np.sin( data["index", 'cylindrical_theta'][idx] ) * vrot
      return data.ds.arr( vx, "km/s" ).in_base( data.ds.unit_system.name )
   ds.add_field( ("gas", "local_rotational_velocity_x"), function=_local_rotational_velocity_x,
                 sampling_type="cell", take_log=False, particle_type=False, units="km/s" )

   def _local_rotational_velocity_y( field, data ):
      vy = np.zeros( data[("gas", "velocity_y")].shape )
      for r, vrot in zip(radius, gas_vrot):
         idx = np.where( (data[("index", "cylindrical_r")].in_units("kpc") >= (r - 0.5*dr)) &
                         (data[("index", "cylindrical_r")].in_units("kpc") <  (r + 0.5*dr)) )
         vy[idx] =  np.cos( data["index", 'cylindrical_theta'][idx] ) * vrot
      return data.ds.arr( vy, "km/s" ).in_base( data.ds.unit_system.name )
   ds.add_field( ("gas", "local_rotational_velocity_y"), function=_local_rotational_velocity_y,
                 sampling_type="cell", take_log=False, particle_type=False, units="km/s" )

   def _velocity_minus_local_rotational_velocity_squared( field, data ):
           return ( data[("gas", "velocity_x")] - data[("gas", "local_rotational_velocity_x")] )**2 + \
                  ( data[("gas", "velocity_y")] - data[("gas", "local_rotational_velocity_y")] )**2 + \
                  ( data[("gas", "velocity_z")]                                                )**2
   ds.add_field( ("gas", "velocity_minus_local_rotational_velocity_squared"), function=_velocity_minus_local_rotational_velocity_squared,
                 sampling_type="cell", take_log=False, particle_type=False, units="km**2/s**2" )

   prof     = yt.ProfilePlot( sp_dense, ("index", "cylindrical_r"),  ("gas", "velocity_minus_local_rotational_velocity_squared"),
                              weight_field=("gas", "cell_mass"), n_bins=nbin, x_log=False )
   gas_vdis = np.sqrt( prof.profiles[0]["velocity_minus_local_rotational_velocity_squared"] ).in_units("km/s").d



#  plot
#  ==================================================================================
   f, ax = plt.subplots( 2, 2 )
#  f.subplots_adjust( hspace=0.0, wspace=0.0 )
   f.subplots_adjust( wspace=0.4 )

   [ f.axes[t].set_xlim( 0.0, 14.0 ) for t in range(0,4,1) ]

   ax[1][0].set_xlabel( "$\mathrm{Cylindrical\ radius\ [kpc]}$", fontsize='large' )
   ax[1][1].set_xlabel( "$\mathrm{Cylindrical\ radius\ [kpc]}$", fontsize='large' )

#  (1) gas surface density
   ax[0][0].plot( radius, gas_dens, 'r-o', lw=2, mec='none', ms=markersize )
   ax[0][0].set_yscale( 'log', nonposy='clip' )
   ax[0][0].set_ylim( 1.0e0, 2.0e3 )
   ax[0][0].yaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[2.0,5.0,8.0]) )
   ax[0][0].set_ylabel( "$\mathrm{\Sigma_{gas}\ [M_{\odot}/pc^2]}$", fontsize='large' )

#  (2) gas temperature
   ax[0][1].plot( radius, gas_temp, 'r-o', lw=2, mec='none', ms=markersize )
   ax[0][1].set_yscale( 'log', nonposy='clip' )
   ax[0][1].set_ylim( 1.0e1, 1.0e4 )
   ax[0][1].yaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[2.0,5.0,8.0]) )
   ax[0][1].set_ylabel( "$\mathrm{T\ [K]}$", fontsize='large' )

#  (3) gas rotational velocity
   ax[1][0].plot( radius, gas_vrot, 'r-o', lw=2, mec='none', ms=markersize )
   ax[1][0].set_ylim( 50.0, 280.0 )
   ax[1][0].set_ylabel( "$\mathrm{v_{rot,gas}\ [km/s]}$", fontsize='large' )

#  (4) gas velocity dispersion
   ax[1][1].plot( radius, gas_vdis, 'r-o', lw=2, mec='none', ms=markersize )
   ax[1][1].set_ylim( 00.0, 100.0 )
   ax[1][1].set_ylabel( "$\mathrm{\sigma_{gas}\ [km/s]}$", fontsize='large' )

#  add title
   time = ds.current_time.in_units('Myr')
   plt.suptitle( "t = %6.2f %s"%(time.d, time.units), fontsize='large' )

#  show/save figure
   plt.savefig( fileout+'_'+ds.basename+".png", bbox_inches='tight', pad_inches=0.05, dpi=dpi )
#  plt.show()

