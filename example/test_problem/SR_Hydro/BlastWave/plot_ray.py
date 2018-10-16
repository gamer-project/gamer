import argparse
import sys
import yt
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

####################  ON-DISK DATA  ###############################
# define primitive density field
def _density_pri( field, data):
   return data["Dens_Pri"]*ds.mass_unit/ds.length_unit**3

# define 4-velocity in x-direction
def _4velocity_x ( field, data):
   return data["Ux"]*ds.length_unit/ds.time_unit

# define 4-velocity in y-direction
def _4velocity_y ( field, data):
   return data["Uy"]*ds.length_unit/ds.time_unit

# define 4-velocity in z-direction
def _4velocity_z ( field, data):
   return data["Uz"]*ds.length_unit/ds.time_unit

# define pressure field
def _pressure_sr( field, data ):
   return data["Pres"]*ds.mass_unit/(ds.length_unit*ds.time_unit**2)

####################  DERIVED DATA  ############################

# define Lorentz factor
def _lorentz_factor (field, data):
   return np.sqrt(1+data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2)

# define magnitude of 4-velocity
def _4velocity_mag ( field, data):
   return np.sqrt( data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2 )*ds.length_unit/ds.time_unit

# define 3-velocity in x-direction
def _3velocity_x ( field, data):
   factor = 1/np.sqrt(1 + data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2 )
   return factor * data["Ux"]*ds.length_unit/ds.time_unit

# define 3-velocity in y-direction
def _3velocity_y ( field, data):
   factor = 1/np.sqrt(1 + data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2 )
   return factor * data["Uy"]*ds.length_unit/ds.time_unit

# define 3-velocity in z-direction
def _3velocity_z ( field, data):
   factor = 1/np.sqrt(1 + data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2 )
   return factor * data["Uz"]*ds.length_unit/ds.time_unit

# define magnitude of 3-velocity field
def _3velocity_mag ( field, data):
   factor = 1/np.sqrt(1 + data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2 )
   Vx = factor * data["Ux"]
   Vy = factor * data["Uy"]
   Vz = factor * data["Uz"]
   return np.sqrt(Vx**2 + Vy**2 + Vz**2)*ds.length_unit/ds.time_unit

# define conserved density field
def _density_cons( field, data):
   factor = np.sqrt(1 + data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2 )
   return data["Dens_Pri"]*factor*ds.mass_unit/ds.length_unit**3

## define momentum in x-direction
#def _momentum_x_sr ( field, data):
#   factor_1 = np.sqrt(1 + data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2 )
#   factor_2 = data.ds.gamma / ( data.ds.gamma - 1.0 )
#   factor_3 = data["Dens_Pri"] + factor_2 * data["Pres"]
#   factor_4 = factor_3 * factor_1
#   return factor_4 * data["Ux"]*ds.mass_unit/(ds.time_unit*ds.length_unit**2)
#
## define momentum in y-direction
#def _momentum_y_sr ( field, data):
#   factor_1 = np.sqrt( 1 + data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2 )
#   factor_2 = data.ds.gamma / ( data.ds.gamma - 1.0 )
#   factor_3 = data["Dens_Pri"] + factor_2 * data["Pres"]
#   factor_4 = factor_3 * factor_1
#   return factor_4 * data["Uy"]*ds.mass_unit/(ds.time_unit*ds.length_unit**2)
#
## define momentum in z-direction
#def _momentum_z_sr ( field, data):
#   factor_1 = np.sqrt( 1 + data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2 )
#   factor_2 = data.ds.gamma / ( data.ds.gamma - 1.0 )
#   factor_3 = data["Dens_Pri"] + factor_2 * data["Pres"]
#   factor_4 = factor_3 * factor_1
#   return factor_4 * data["Uz"]*ds.mass_unit/(ds.time_unit*ds.length_unit**2)
#
## define conserved energy field
#def _energy_cons( field, data):
#   factor_0 = 1 + data["Ux"]**2 + data["Uy"]**2 + data["Uz"]**2
#   factor_2 = data.ds.gamma / ( data.ds.gamma - 1.0 )
#   factor_3 = data["Dens_Pri"] + factor_2 * data["Pres"]
#   return (factor_3 * factor_0 - data["Pres"] )*(ds.mass_unit/(ds.length_unit*ds.time_unit**2))


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot sr-hydro velocity slices for the blast wave test' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False,  type=str, dest='prefix',
                     help='data path prefix [%(default)s]', default='./' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]) ),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap    = 'arbre'

field       = 'density'    # to change the target field, one must modify set_unit() accordingly
#field       = 'Ux'    # to change the target field, one must modify set_unit() accordingly
#field       = 'Uy'    # to change the target field, one must modify set_unit() accordingly
#field       = 'Uz'    # to change the target field, one must modify set_unit() accordingly
#field       = 'Pres'           # to change the target field, one must modify set_unit() accordingly
#field       = '4velocity_mag'  # to change the target field, one must modify set_unit() accordingly
#field       = '3velocity_x'    # to change the target field, one must modify set_unit() accordingly
#field       = '3velocity_y'    # to change the target field, one must modify set_unit() accordingly
#field       = '3velocity_z'    # to change the target field, one must modify set_unit() accordingly
#field       = '3velocity_mag'  # to change the target field, one must modify set_unit() accordingly
#field       = 'density_cons'   # to change the target field, one must modify set_unit() accordingly
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

# add new derived field
   ds.add_field( ("gamer", "Lorentz_factor") , function=_density_pri  , sampling_type="cell", units="" )
   ds.add_field( ("gamer", "density")      , function=_density_pri  , sampling_type="cell", units="code_mass/code_length**3" )
   ds.add_field( ("gamer", "4-velocity_x") , function=_4velocity_x  , sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "4-velocity_y") , function=_4velocity_y  , sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "4-velocity_z") , function=_4velocity_z  , sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "pressure")     , function=_pressure_sr  , sampling_type="cell", units="code_mass/(code_length*code_time**2)" )
   ds.add_field( ("gamer", "4-velocity_z") , function=_4velocity_z  , sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "magnitude of 4-velocity"), function=_4velocity_mag, sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "3-velocity_x")  , function=_3velocity_x  , sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "3-velocity_y")  , function=_3velocity_y  , sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "3-velocity_z")  , function=_3velocity_z  , sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "magnitude of 3-velocity"), function=_3velocity_mag, sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "conservative density") , function=_density_cons , sampling_type="cell", units="code_mass/code_length**3" )
#   ds.add_field( ("gamer", "momentum_x_sr"), function=_momentum_x_sr, sampling_type="cell", units="code_mass/(code_time*code_length**2)" )
#   ds.add_field( ("gamer", "momentum_y_sr"), function=_momentum_y_sr, sampling_type="cell", units="code_mass/(code_time*code_length**2)" )
#   ds.add_field( ("gamer", "momentum_z_sr"), function=_momentum_z_sr, sampling_type="cell", units="code_mass/(code_time*code_length**2)" )
#   ds.add_field( ("gamer", "energy_cons")  , function=_energy_cons  , sampling_type="cell", units="code_mass/(code_length*code_time**2)" )

   my_ray = ds.ray((-0, 0.5, 0.5), (1.0, 0.5, 0.5))
   srt = np.argsort(my_ray['x'])
   density = my_ray[field][srt]
   
#   plt.semilogy(np.array(my_ray['x'][srt]), np.array(my_ray[field][srt]), 'bo')
   plt.plot(np.array(my_ray['x'][srt]), np.array(my_ray[field][srt]), 'bo', markersize=2)
   plt.xlabel('x')
   plt.ylabel(field)
   plt.savefig("density_xsweep.png")
