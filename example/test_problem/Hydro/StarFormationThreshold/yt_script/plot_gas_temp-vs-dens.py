import argparse
import sys
import yt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas density-temperature phase diagram' )

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
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = '../'

colormap    = 'algae'
nbin        = 128
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

x_lim_min = 1.0e-29
x_lim_max = 1.0e-21
y_lim_min = 1.0e0
y_lim_max = 1.0e8
z_lim_min = 1.0e0
z_lim_max = 1.0e8

for ds in ts.piter():

#  plot
   temp_dens = yt.PhasePlot( ds, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'cell_mass'),
                             weight_field=None, x_bins=nbin, y_bins=nbin )

   temp_dens.set_unit( 'cell_mass', 'Msun' )
   temp_dens.set_xlim( x_lim_min, x_lim_max )
   temp_dens.set_ylim( y_lim_min, y_lim_max )
   temp_dens.set_zlim( ('gas', 'cell_mass'), z_lim_min, z_lim_max )
   temp_dens.set_cmap( ('gas', 'cell_mass'), colormap )
   temp_dens.set_colorbar_label( ('gas', 'cell_mass'), "Mass ($\mathrm{M}_{\odot}$)" )
   temp_dens.annotate_text( xpos=x_lim_min*10**(0.80*(np.log10(x_lim_max)-np.log10(x_lim_min))),
                            ypos=y_lim_min*10**(0.95*(np.log10(y_lim_max)-np.log10(y_lim_min))),
                            text='$t$ = {:.1f} {:s}'.format( ds.current_time.in_units('Myr').d, 'Myr' ),
                            color='white' )
   temp_dens.save( mpl_kwargs={'dpi':dpi} )
