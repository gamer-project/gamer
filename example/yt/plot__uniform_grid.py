
# ref: https://yt-project.org/docs/dev/examining/generic_array_data.html

import yt
import sys
import numpy as np


# input parameters
file_in     = 'filename_in'            # input filename
file_out    = 'filename_out'           # output filename
field       = 'field_name'             # target field
dim         = 2                        # data dimensionality (2/3)
N           = 512                      # data dimensions = N^dim
float_type  = 'float32'                # data floating-point type (float32/64)
dh          = 1.23456789               # cell width in length_unit
length_unit = (1.0, 'kpc')             # unit of dh
data_unit   = (1.0, 'Msun/kpc**2')     # unit of input data
plot_unit   = 'Msun/kpc**2'            # unit of output image
dpi         = 150                      # image resolution in dot per inch


# load data
arr1d = np.fromfile( file_in, dtype=float_type )

if   dim == 2:
   arr3d = np.transpose( arr1d.reshape((1,N,N)) )
   bbox  = np.array( [ [0,N*dh], [0,N*dh], [0,dh] ] )
elif dim == 3:
   arr3d = np.transpose( arr1d.reshape((N,N,N)) )
   bbox = np.array( [ [0,N*dh], [0,N*dh], [0,N*dh] ] )
else:
   sys.exit( 'ERROR : dim must be 2/3 !!' )

data = {}
data[field] = ( arr3d, str(data_unit[0])+"*"+data_unit[1] )

ds = yt.load_uniform_grid( data, arr3d.shape, length_unit=length_unit, bbox=bbox, nprocs=1 )


# plot
sz = yt.SlicePlot( ds, "z", field, center='c', origin=('left','window') )
#sz = yt.ProjectionPlot( ds, "z", field, center='c', origin=('left','window') )
sz.set_unit( field, plot_unit )

#sz.set_zlim( field, 1.0e4, 1.0e8 )
#sz.set_colorbar_label( field, "Surface density [$\\mathrm{M_\\odot/kpc^2}$]" )
#sz.annotate_scale( coeff=5.0, unit='kpc', pos=(0.15,0.10) )
#sz.hide_axes()
#sz.set_cmap( field, "afmhot" )
#sz.set_figure_size( 16 )
#sz.set_buff_size( 2048 )
#sz.set_font( {'size':36} )

sz.save( file_out+".png", mpl_kwargs={"dpi":dpi} )

