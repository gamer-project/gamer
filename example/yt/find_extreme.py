import yt
import sys


# parameters
field = 'Dens'
file  = 'Data_000000'
mode  = 'max'           # max/min

# load data
ds = yt.load( file )

# show the field list
#ds.field_list

# fine the max/min and location of the target field
if   mode == 'max':
   val, center = ds.find_max( field )
elif mode == 'min':
   val, center = ds.find_min( field )
else:
   sys.exit( 'ERROR : mode must be max/min !!' )
