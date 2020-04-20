
# ref: https://yt-project.org/docs/dev/examining/generic_array_data.html

import yt
import h5py
import numpy as np


# user-specified parameters
file_in = 'Cube_x0.000-3.000_y0.000-3.000_z0.000-3.000_lv1.hdf5'  # input filename
fields  = [ 'Dens', ]                                             # target field(s)


# load data
f         = h5py.File( file_in, 'r' )
mode      = f['Info']['OutputMode']
dimension = f['Info']['GridDimension']
time      = f['Info']['Time']
dh        = f['Info']['CellWidth']
box_size  = f['Info']['SubdomainSize']
left_edge = f['Info']['SubdomainLeftEdge']
unit_l    = ( f['Info']['Unit_L'], 'cm' )
unit_m    = ( f['Info']['Unit_M'], 'g' )
unit_t    = ( f['Info']['Unit_T'], 's' )

# overwrite the subdomain coordinates for the projected data
if mode >= 4 and mode <=6:
   dir            = mode - 4
   left_edge[dir] = left_edge[dir] + 0.5*( box_size[dir] - dh )
   box_size [dir] = dh

bbox  = np.array( [ [left_edge[0], left_edge[0]+box_size[0]],
                    [left_edge[1], left_edge[1]+box_size[1]],
                    [left_edge[2], left_edge[2]+box_size[2]] ] )
units = [ f['Data'][k].attrs['Unit'].decode('utf-8') for k in fields ]
data  = { k:(f['Data'][k][()].transpose(),u) for k,u in zip(fields,units) }

ds = yt.load_uniform_grid( data=data, domain_dimensions=dimension,
                           length_unit=unit_l, mass_unit=unit_m, time_unit=unit_t,
                           bbox=bbox, nprocs=1, sim_time=time,
                           periodicity=(False,False,False) )


# plot
sz = yt.SlicePlot( ds, 'z', fields, center='c' )
sz.save()
