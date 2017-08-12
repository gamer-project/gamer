import yt

yt.enable_parallelism()

idx_start   = 0
idx_end     = 50
didx        = 1
prefix      = '../'

colormap    = 'algae'
width_kpc   = 30
nbin        = 300
dpi         = 150

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
#ts = yt.load( 'Data_??????' )

for ds in ts.piter():

#  define center as the location of peak gas density within 1 kpc from the center of gas mass
   v, cen1 = ds.h.find_max( ("gas", "density") )
   sp1  = ds.sphere( cen1, (30.0, "kpc") )
   cen2 = sp1.quantities.center_of_mass( use_gas=True, use_particles=False ).in_units( "kpc" )
   sp2  = ds.sphere( cen2, (1.0, "kpc") )
   cen3 = sp2.quantities.max_location( ("gas", "density") )
   cen  = ds.arr( [cen3[1].d, cen3[2].d, cen3[3].d], 'code_length' )


#  only include the data within a sphere with a radius of width_kpc
   sp = ds.sphere( cen, (0.5*width_kpc, "kpc") )


#  plot
   temp_dens = yt.PhasePlot( sp, ("gas", "density"), ("gas", "temperature"), ("gas", "cell_mass"),
                             weight_field=None, x_bins=nbin, y_bins=nbin )
   temp_dens.set_unit( "cell_mass", "Msun" )
   temp_dens.set_xlim( 1.0e-29, 1.0e-21 )
   temp_dens.set_ylim( 1.0e1, 1.0e7 )
   temp_dens.set_zlim( ("gas", "cell_mass"), 1.0e3, 1.0e8 )
   temp_dens.set_cmap( ("gas", "cell_mass"), colormap )
   temp_dens.set_colorbar_label( ("gas", "cell_mass"), "Mass ($\mathrm{M}_{\odot}$)" )
#  temp_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   temp_dens.save( mpl_kwargs={"dpi":dpi} )


