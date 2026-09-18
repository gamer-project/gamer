import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yt
from mpl_toolkits.axes_grid1 import AxesGrid
import argparse
import sys

#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Get slices' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )

field       = ['velocity',      'vorticity',   'magnetic_field' ]
field_unit  = ['code_velocity', '1/code_time', 'code_magnetic'  ]
colormap    = ['viridis',       'cividis',     'plasma'         ]
zmin        = [ 1.0e-4,          1.0e-1,        1.0e-4          ]
zmax        = [ 1.0,             1.0e+2,        1.0             ]


dpi       = 150
fontsize  = 24
titlepad  = 16

yt.enable_parallelism()
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
   nfield   = 3   if ('gas', 'magnetic_field_magnitude') in ds.derived_field_list else 2
   figwidth = 3.2 if ('gas', 'magnetic_field_magnitude') in ds.derived_field_list else 2.2

#  plot
   fig = plt.figure()
   fig.dpi = dpi
   grid = AxesGrid( fig, (0.1, 0.05, figwidth, 2.7), nrows_ncols=(1, nfield), axes_pad=(3,0.5), label_mode="all",
                          share_all=True, cbar_location="right", cbar_mode="each", cbar_size="2%", cbar_pad="2%")

   slc = [None]*nfield
   for i in range(nfield):
      fieldname = field[i]+"_magnitude"
      slc[i] = yt.SlicePlot( ds, 0, fields = fieldname, center = 'c')
      slc[i].set_background_color( fieldname )
      slc[i].set_axes_unit( 'code_length' )
      slc[i].set_unit( fieldname, field_unit[i])
      slc[i].set_zlim( fieldname, zmin[i], zmax[i] )
      slc[i].set_cmap( fieldname, colormap[i] )
      slc[i].set_font( {'size':fontsize} )
      slc[i].annotate_grids( periodic=False )
#     slc[i].annotate_streamlines(("gas", field[i]+"_x"), ("gas", fiel2[i]+"_y"), color='black')
      slc[i].annotate_quiver(("gas", field[i]+"_x"), ("gas", field[i]+"_y"), color='black')

      plot = slc[i].plots[fieldname]
      plot.figure = fig
      plot.axes = grid[i].axes
      plot.cax = grid.cbar_axes[i]
      slc[i]._setup_plots()
      grid[i].set_title(fieldname, fontsize=fontsize, pad=titlepad)

   fig.savefig("fig_%s_Slice.png"%(ds), bbox_inches='tight',pad_inches=0.02)


