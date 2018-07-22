import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot gravitational potential errors on a slice ' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='filename_in',
                     help='input filename [%(default)s]', default='PotError.bin' )

args        = parser.parse_args()
filename_in = args.filename_in


# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]) ),
print( '' )
print( '-------------------------------------------------------------------\n' )


# input parameters
colormap     = 'arbre'
field        = 'MomY'   # gas fields have been overwritten: [MomX/Y]=absolute/relative errors in potential
center_mode  = 'c'
width        = 2.0e1
dpi          = 150
filename_out = 'Fig__PotError_Slice.png'


# plot
ds = yt.load( filename_in )
sz = yt.SlicePlot( ds, 'z', field, center=center_mode  )
sz.set_width( width, 'code_length' )
#sz.set_zlim( field, 1.0e-3, 1.0e-2 )
sz.set_log( field, True )
sz.set_cmap( field, colormap )
sz.set_axes_unit( 'code_length' )
sz.set_colorbar_label( field, 'Relative potential error' )
sz.annotate_grids( periodic=False )
sz.save( filename_out, mpl_kwargs={"dpi":dpi} )

