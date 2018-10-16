import argparse
import sys
import yt
import math

# define Q field 
def _Q( field, data ):
   return data["Q"]*ds.mass_unit/ds.length_unit**3

# define magnitude of 4-velocity
def _4velocity_mag( field, data ):
   if data["Q"] == data["Q"]:
     Ux = data["MomX"]/data["Q"]
     Uy = data["MomY"]/data["Q"]
     Uz = data["Momz"]/data["Q"]
     return math.sqrt(Ux**2 + Uy**2 + Uz**2)
   else:
     return 0
   

# define 4-velocity in x-direction
def _4velocity_x( field, data ):
   if data["Q"] == data["Q"]:
     Ux = data["MomX"]/data["Q"]
     return Ux
   else:
     return 0

# define 4-velocity in y-direction
def _4velocity_y( field, data ):
   if data["Q"] == data["Q"]:
     Uy = data["MomY"]/data["Q"]
     return Uy
   else:
     return 0

# define 4-velocity in z-direction
def _4velocity_z( field, data ):
   if data["Q"] == data["Q"]:
     Uz = data["MomZ"]/data["Q"]
     return Uz
   else:
     return 0

# define magnitude of 3-velocity
def _3velocity_mag( field, data ):
   if data["Q"] == data["Q"]:
     Ux = data["MomX"]/data["Q"]
     Uy = data["MomY"]/data["Q"]
     Uz = data["Momz"]/data["Q"]
     c = 1*ds.length_unit/ds_time_unit
     factor = math.sqrt(1+(Ux/c)**2+(Uy/c)**2+(Uz/c)**2)
     Vx = Ux / factor
     Vy = Uy / factor
     Vz = Uz / factor
     return math.sqrt(Vx**2 + Vy**2 + Vz**2)
   else:
     return 0

# define 3-velocity in x-direction
def _3velocity_x( field, data ):
   if data["Q"] == data["Q"]:
     Ux = data["MomX"]/data["Q"]
     Uy = data["MomY"]/data["Q"]
     Uz = data["Momz"]/data["Q"]
     c = 1*ds.length_unit/ds_time_unit
     factor = math.sqrt(1+(Ux/c)**2+(Uy/c)**2+(Uz/c)**2)
     return Ux / factor
   else:
     return 0

# define 3-velocity in y-direction
def _3velocity_y( field, data ):
   if data["Q"] == data["Q"]:
     Ux = data["MomX"]/data["Q"]
     Uy = data["MomY"]/data["Q"]
     Uz = data["Momz"]/data["Q"]
     c = 1*ds.length_unit/ds_time_unit
     factor = math.sqrt(1+(Ux/c)**2+(Uy/c)**2+(Uz/c)**2)
     return Uy / factor
   else:
     return 0

# define 3-velocity in z-direction
def _3velocity_z( field, data ):
   if data["Q"] == data["Q"]:
     Ux = data["MomX"]/data["Q"]
     Uy = data["MomY"]/data["Q"]
     Uz = data["Momz"]/data["Q"]
     c = 1*ds.length_unit/ds_time_unit
     factor = math.sqrt(1+(Ux/c)**2+(Uy/c)**2+(Uz/c)**2)
     return Uz / factor
   else:
     return 0

# define primitive density
def _density_pri( field, data ):
   if data["Q"] == data["Q"]:
     Ux = data["MomX"]/data["Q"]
     Uy = data["MomY"]/data["Q"]
     Uz = data["Momz"]/data["Q"]
     c = 1*ds.length_unit/ds_time_unit
     factor = math.sqrt(1+(Ux/c)**2+(Uy/c)**2+(Uz/c)**2)
     return data["Dens"] / factor
   else:
     return data["Dens"]

# define pressure
def _pressure( field, data ):
   if data["Q"] == data["Q"]:
     Ux = data["MomX"]/data["Q"]
     Uy = data["MomY"]/data["Q"]
     Uz = data["Momz"]/data["Q"]
     c = 1*ds.length_unit/ds_time_unit
     factor = math.sqrt(1+(Ux/c)**2+(Uy/c)**2+(Uz/c)**2)
     pres = (ds.gamma-1/ds.gamma)*( (data["Q"] - data["Dens"]) / factor )*c**2
     return pres
   else:
     pres = (ds.gamma-1)*(data["Engy"] - data["Dens"]*c**2)
     return pres







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
field       = 'pressure_pri'    # to change the target field, one must modify set_unit() accordingly
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

# add new derived field
   ds.add_field( ("gamer", "4velocity_mag"), function=_4velocity_mag, sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "4velocity_x"),   function=_4velocity_x,   sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "4velocity_y"),   function=_4velocity_y,   sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "4velocity_z"),   function=_4velocity_z,   sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "4velocity_mag"), function=_3velocity_mag, sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "3velocity_x"),   function=_3velocity_x,   sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "3velocity_y"),   function=_3velocity_y,   sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "3velocity_z"),   function=_3velocity_z,   sampling_type="cell", units="code_length/code_time" )
   ds.add_field( ("gamer", "density_pri"),   function=_density_pri,   sampling_type="cell", units="code_mass/code_length**3" )
   ds.add_field( ("gamer", "pressure_pri"),  function=_pressure,      sampling_type="cell", units="code_mass/(code_length*code_time**2)" )

   sz = yt.SlicePlot( ds, 'z', field, center_mode  )
   sz.set_zlim( field, -1.0, 1.0 )
   sz.set_log( field, False )
   sz.set_cmap( field, colormap )
#   sz.set_unit( field, 'code_mass/(code_time*code_length**2)' )
   sz.set_unit( field, 'code_mass/(code_length*code_time**2)' )
   sz.set_axes_unit( 'code_length' )
   sz.annotate_timestamp( time_unit='code_time', corner='upper_right', time_format='t = {time:.4f} {units}' )
   sz.annotate_grids( periodic=False )
   sz.save( mpl_kwargs={"dpi":dpi} )
