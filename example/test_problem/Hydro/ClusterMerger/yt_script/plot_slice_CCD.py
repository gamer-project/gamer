#====================================================================================================
# Import
#====================================================================================================
import argparse
import sys
import yt
from yt.units import Mpc
import h5py

sys.dont_write_bytecode = True
import MergerUtilities as mus



#====================================================================================================
# Consant
#====================================================================================================
# some constants in CGS
amu  = 1.660539040e-24
mu   = 0.59242761692650336
ne   = 0.52                # electron number fraction
keV  = 1.6021766208e-9
kB   = 1.38064852e-16
msun = 1.9885e33
kpc  = 3.08567758149e21



#====================================================================================================
# Main
#====================================================================================================
# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the density and entropy slices of CCD simulations' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
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
print( ' '.join(sys.argv) )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix


def _entropy(field, data):
    return data["gas", "temperature"]*kB/keV*(data["gas", "density"]*ne/(amu*mu))**(-2.0/3.0)
yt.add_field(("gas", "entropy"),function=_entropy,sampling_type="local",units="K*cm**2/g**(2/3)")

# define the particle filters for SMBH particles
def smbh1( pfilter, data ):
   filter = (data[ "all", "ParType" ] == 4) & (data["all", "ParHalo"] == 0)
   return filter

def smbh2( pfilter, data ):
   filter = (data[ "all", "ParType" ] == 4) & (data["all", "ParHalo"] == 1)
   return filter

yt.add_particle_filter( "smbh1", function=smbh1, filtered_type="all", requires=["ParType"] )
yt.add_particle_filter( "smbh2", function=smbh2, filtered_type="all", requires=["ParType"] )


yt.enable_parallelism()
ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )


# Get clusters center
cenX0, cenY0, cenZ0 = mus.getClusterCen(prefix, 0, idx_start, idx_end, didx)
cenX1, cenY1, cenZ1 = mus.getClusterCen(prefix, 1, idx_start, idx_end, didx)


colormap_dens = 'algae'
colormap_temp = 'afmhot'
colormap_cdm  = 'arbre'
center_mode   = 'c'
dpi           = 150

plot_width_one = (200,  'kpc')
plot_width_all = (5000, 'kpc')


N0 = len(cenX0)
N1 = len(cenX1)
ii = -1
for ds in ts.piter():
   ii += 1
   # apply the particle filters to the dataset
   ds.add_particle_filter( 'smbh1' )
   ds.add_particle_filter( 'smbh2' )

   # all
   field   = 'density'
   sz_dens = yt.SlicePlot( ds, 'z', field, center=center_mode, width=plot_width_all )

   sz_dens.set_zlim( field, 1.0e-27, 1.0e-22 )
   sz_dens.set_cmap( field, colormap_dens )
   sz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_dens.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='b', ptype='smbh1' )
   sz_dens.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='r', ptype='smbh2' )
   sz_dens.save( prefix+"Data_%06d_density_all.png"%(ii), mpl_kwargs={"dpi":dpi} )

   field   = 'entropy'
   sz_entr = yt.SlicePlot( ds, 'z', field, center=center_mode, width=plot_width_all )
   sz_entr.set_zlim( field, 1.0e+1, 1.0e+5 )
   sz_entr.set_cmap( field, colormap_dens )
   sz_entr.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_entr.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='b', ptype='smbh1' )
   sz_entr.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='r', ptype='smbh2' )
   sz_entr.save( prefix+"Data_%06d_entropy_all.png"%(ii), mpl_kwargs={"dpi":dpi} )

   # major
   field   = 'density'
   sz_dens = yt.SlicePlot( ds, 'z', field, center=[cenX0[ii], cenY0[ii], cenZ0[ii]], width=plot_width_one )

   sz_dens.set_zlim( field, 1.0e-27, 1.0e-22 )
   sz_dens.set_cmap( field, colormap_dens )
   sz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_dens.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='b', ptype='smbh1' )
   sz_dens.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='r', ptype='smbh2' )
   sz_dens.save( prefix+"Data_%06d_density_major.png"%(ii), mpl_kwargs={"dpi":dpi} )

   field   = 'entropy'
   sz_entr = yt.SlicePlot( ds, 'z', field, center=[cenX0[ii], cenY0[ii], cenZ0[ii]], width=plot_width_one )
   sz_entr.set_zlim( field, 1.0e+1, 1.0e+5 )
   sz_entr.set_cmap( field, colormap_dens )
   sz_entr.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_entr.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='b', ptype='smbh1' )
   sz_entr.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='r', ptype='smbh2' )
   sz_entr.save( prefix+"Data_%06d_entropy_major.png"%(ii), mpl_kwargs={"dpi":dpi} )

   if ii >= N1: continue

   # minor
   field   = 'density'
   sz_dens = yt.SlicePlot( ds, 'z', field, center=[cenX1[ii], cenY1[ii], cenZ1[ii]], width=plot_width_one )

   sz_dens.set_zlim( field, 1.0e-27, 1.0e-22 )
   sz_dens.set_cmap( field, colormap_dens )
   sz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_dens.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='b', ptype='smbh1' )
   sz_dens.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='r', ptype='smbh2' )
   sz_dens.save( prefix+"Data_%06d_density_minor.png"%(ii), mpl_kwargs={"dpi":dpi} )

   field   = 'entropy'
   sz_entr = yt.SlicePlot( ds, 'z', field, center=[cenX1[ii], cenY1[ii], cenZ1[ii]], width=plot_width_one )
   sz_entr.set_zlim( field, 1.0e+1, 1.0e+5 )
   sz_entr.set_cmap( field, colormap_dens )
   sz_entr.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_entr.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='b', ptype='smbh1' )
   sz_entr.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='r', ptype='smbh2' )
   sz_entr.save( prefix+"Data_%06d_entropy_minor.png"%(ii), mpl_kwargs={"dpi":dpi} )


