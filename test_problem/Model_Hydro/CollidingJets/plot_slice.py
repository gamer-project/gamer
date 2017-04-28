import yt

yt.enable_parallelism()

idx_start = 0
idx_end   = 50
didx      = 1
prefix    = 'Data_'

ts = yt.load( [ prefix+'%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
#ts = yt.load( 'Data_??????' )

for ds in ts.piter():

   sz = yt.SlicePlot( ds, 'z', 'temperature', center='c', width=[40,10] )
   sz.set_unit( 'temperature', 'keV', equivalency='thermal' )
   sz.set_zlim( 'temperature', 1.0e0, 2.0e2 )
   sz.set_cmap( 'temperature', 'afmhot' )
   sz.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz.save( mpl_kwargs={"dpi":300} )
