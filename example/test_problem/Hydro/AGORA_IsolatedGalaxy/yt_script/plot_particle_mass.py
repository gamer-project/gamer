import yt

yt.enable_parallelism()

idx_start   = 0
idx_end     = 50
didx        = 1
prefix      = '../Data_'

colormap    = 'algae'
width_kpc   = 30
depth_kpc   = 30
center_mode = 'c'
dpi         = 150

ts = yt.load( [ prefix+'%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
#ts = yt.load( 'Data_??????' )


# define the particle filter for the newly formed stars
def new_star( pfilter, data ):
   filter = data[ "all", "ParCreTime" ] > 0
   return filter

yt.add_particle_filter( "new_star", function=new_star, filtered_type="all", requires=["ParCreTime"] )

AllPar = ( 'all',      'particle_mass' )
NewPar = ( 'new_star', 'particle_mass' )


for ds in ts.piter():

#  add the particle filter
   ds.add_particle_filter( "new_star" )


#  face-on (all particles)
   pz = yt.ParticleProjectionPlot( ds, 'z', AllPar, center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   pz.set_unit( AllPar, 'Msun' )
   pz.set_zlim( AllPar, 1.0e6, 1.0e8 )
   pz.set_cmap( AllPar, colormap )
   pz.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   pz.save( name=ds.basename+'_AllPar', mpl_kwargs={"dpi":dpi} )


#  face-on (new particles)
   pz = yt.ParticleProjectionPlot( ds, 'z', NewPar, center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   pz.set_unit( NewPar, 'Msun' )
#  pz.set_zlim( NewPar, 1.0e6, 1.0e8 )
   pz.set_cmap( NewPar, colormap )
   pz.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   pz.save( name=ds.basename+'_NewPar', mpl_kwargs={"dpi":dpi} )


#  edge-on (all particles)
   px = yt.ParticleProjectionPlot( ds, 'x', AllPar, center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   px.set_unit( AllPar, 'Msun' )
   px.set_zlim( AllPar, 1.0e6, 1.0e8 )
   px.set_cmap( AllPar, colormap )
   px.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   px.save( name=ds.basename+'_AllPar', mpl_kwargs={"dpi":dpi} )


#  edge-on (new particles)
   px = yt.ParticleProjectionPlot( ds, 'x', NewPar, center=center_mode, width=(width_kpc,'kpc'), depth=(depth_kpc,'kpc') )
   px.set_unit( NewPar, 'Msun' )
#  px.set_zlim( NewPar, 1.0e6, 1.0e8 )
   px.set_cmap( NewPar, colormap )
   px.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   px.save( name=ds.basename+'_NewPar', mpl_kwargs={"dpi":dpi} )
