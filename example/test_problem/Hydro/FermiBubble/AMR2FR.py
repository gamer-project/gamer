#====================================================================================================
# Imports
#====================================================================================================
import numpy as np
import h5py
import yt
import time



#====================================================================================================
# Constants
#====================================================================================================
PARTICLE_MASS = 1.6726231e-24 # proton mass in gram



#====================================================================================================
# Derived fields
#====================================================================================================
def pressure_sr( field, data ):
    from yt.units import speed_of_light_cgs, boltzmann_constant_cgs
    unit_v = speed_of_light_cgs

    pres  = data["frame_density"] * data["Temp"] * yt.YTQuantity(1.0, "K")
    pres /= PARTICLE_MASS * yt.YTQuantity(1.0, "g")
    pres /= speed_of_light_cgs**2
    pres *= unit_v**2
    pres *= boltzmann_constant_cgs
    return pres



#====================================================================================================
# Main
#====================================================================================================
Data  = 'Data_000035'
Path  = './'+Data
level = 4
Field = [ "Dens", "Temp", "Pres", "CRay", "Passive_0001" ]

ds   = yt.load(Path)
dims = ds.domain_dimensions * ds.refine_by**level

ds.add_field( name=("gamer", "Pres"),     function=_pressure, sampling_type="local", units="dyne/cm**2" )
ds.add_field( name=("gamer", "CRay_new"), function=_cray,     sampling_type="local", units="dyne/cm**2" )


# Note that we open with 'w' (write), which will overwrite existing files!
f = h5py.File( "FRB_%s.h5"%Data, mode="w" )

for i in range(len(Field)):
  cube = ds.covering_grid( level, left_edge=ds.domain_left_edge, dims=dims, fields=Field[i] )
  if Field[i] == "CRay":
      f.create_dataset( "Fields/"+Field[i], shape=cube[("gamer", "CRay_new")].shape, data=cube[("gamer", "CRay_new")].astype('<f4') )
  else:
      f.create_dataset( "Fields/"+Field[i], shape=cube[("gamer", Field[i])].shape, data=cube[("gamer", Field[i])].astype('<f4') )


# Store compound data (Info)
numCellX = dims[0]
numCellY = dims[1]
numCellZ = dims[2]

# Create empty record array with 3 rows
ds_dtype = []
ds_dtype.append( ('level',          np.int32  ) )
ds_dtype.append( ('numCellX',       np.int32  ) )
ds_dtype.append( ('numCellY',       np.int32  ) )
ds_dtype.append( ('numCellZ',       np.int32  ) )
ds_dtype.append( ('BoxSizeX',       np.float32) )
ds_dtype.append( ('BoxSizeY',       np.float32) )
ds_dtype.append( ('BoxSizeZ',       np.float32) )
ds_dtype.append( ('Unit_L',         np.float32) )
ds_dtype.append( ('Unit_M',         np.float32) )
ds_dtype.append( ('Unit_T',         np.float32) )
ds_dtype.append( ('Unit_V',         np.float32) )
ds_dtype.append( ('Unit_D',         np.float32) )
ds_dtype.append( ('Unit_E',         np.float32) )
ds_dtype.append( ('Unit_P',         np.float32) )
ds_dtype.append( ('AMRDataName',    '|S11'    ) )
ds_dtype.append( ('EpochTimeStamp', np.int64  ) )

# load list data to record array by field name
ds_arr = np.recarray( (1,), dtype=ds_dtype )
ds_arr['level']          = np.asarray( level            )
ds_arr['numCellX']       = np.asarray( numCellX         )
ds_arr['numCellY']       = np.asarray( numCellY         )
ds_arr['numCellZ']       = np.asarray( numCellZ         )
ds_arr['BoxSizeX']       = np.asarray( ds["BoxSize"][0] )
ds_arr['BoxSizeY']       = np.asarray( ds["BoxSize"][1] )
ds_arr['BoxSizeZ']       = np.asarray( ds["BoxSize"][2] )
ds_arr['Unit_L']         = np.asarray( ds["Unit_L"]     )
ds_arr['Unit_M']         = np.asarray( ds["Unit_M"]     )
ds_arr['Unit_T']         = np.asarray( ds["Unit_T"]     )
ds_arr['Unit_V']         = np.asarray( ds["Unit_V"]     )
ds_arr['Unit_D']         = np.asarray( ds["Unit_D"]     )
ds_arr['Unit_E']         = np.asarray( ds["Unit_E"]     )
ds_arr['Unit_P']         = np.asarray( ds["Unit_P"]     )
ds_arr['AMRDataName']    = Data
ds_arr['EpochTimeStamp'] = np.asarray( time.time_ns()   )

# load data to dataset my_ds1 using recarray
dset = f.create_dataset( 'Info/Keys', data=ds_arr, maxshape=(None) )

f.close()
