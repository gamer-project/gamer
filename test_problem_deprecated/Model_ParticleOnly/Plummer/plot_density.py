import yt

yt.enable_parallelism()

idx_start = 0
idx_end   = 50
didx      = 1
prefix    = 'Data_'

ts = yt.load( [ prefix+'%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   pz_dens = yt.ProjectionPlot( ds, 'z', 'density' )
   pz_dens.set_zlim( 'density', 1.0e-5, 5.0e-2 )
   pz_dens.set_font( {'size':16} )
   pz_dens.annotate_timestamp( corner='upper_right' )
   pz_dens.save()

   pz_Cloud0 = yt.ProjectionPlot( ds, 'z', 'Cloud0' )
   pz_Cloud0.set_zlim( 'Cloud0', 1.0e-5, 5.0e-2 )
   pz_Cloud0.set_font( {'size':16} )
   pz_Cloud0.annotate_timestamp( corner='upper_right' )
   pz_Cloud0.save()

   pz_Cloud1 = yt.ProjectionPlot( ds, 'z', 'Cloud1' )
   pz_Cloud1.set_zlim( 'Cloud1', 1.0e-5, 5.0e-2 )
   pz_Cloud1.set_font( {'size':16} )
   pz_Cloud1.annotate_timestamp( corner='upper_right' )
   pz_Cloud1.save()
