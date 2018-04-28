import yt

yt.enable_parallelism()

idx_start = 0
idx_end   = 10
didx      = 1
prefix    = './'

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   sx_dens = yt.SlicePlot( ds, 'x', 'density', center='c' )
#  sx_dens.set_zlim( 'density', 1.0e-31, 5.0e-27 )
   sx_dens.annotate_timestamp( corner='upper_right' )
   sx_dens.annotate_grids( periodic=True )
   sx_dens.save()
