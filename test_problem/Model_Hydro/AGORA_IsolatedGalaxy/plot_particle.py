import yt

yt.enable_parallelism()

idx_start   = 0
idx_end     = 0
didx        = 1
prefix      = 'Data_'

colormap    = 'algae'
width_kpc   = 30
depth_kpc   = 30
center_mode = 'c'
dpi         = 150

ts = yt.load( [ prefix+'%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
#ts = yt.load( 'Data_??????' )

for ds in ts.piter():

#  face-on
   pz = yt.ParticleProjectionPlot( ds, 'z', 'particle_mass', center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   pz.set_unit( 'particle_mass', 'Msun' )
   pz.set_zlim( 'particle_mass', 1.0e6, 1.0e8 )
   pz.set_cmap( 'particle_mass', colormap )
   pz.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   pz.save( mpl_kwargs={"dpi":dpi} )


#  edge-on
   px = yt.ParticleProjectionPlot( ds, 'x', 'particle_mass', center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   px.set_unit( 'particle_mass', 'Msun' )
   px.set_zlim( 'particle_mass', 1.0e6, 1.0e8 )
   px.set_cmap( 'particle_mass', colormap )
   px.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   px.save( mpl_kwargs={"dpi":dpi} )
