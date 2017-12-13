import yt

yt.enable_parallelism()

idx_start = 0
idx_end   = 10
didx      = 1
prefix    = 'Data_'

ts = yt.load( [ prefix+'%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
#ts = yt.load( 'Data_??????' )

for ds in ts.piter():

   sz_dens = yt.SlicePlot( ds, 'z', 'Dens', center='c' )
#  sz_dens.set_zlim( 'Dens', 1.0e-1, 5.0e5 )
   sz_dens.annotate_timestamp( corner='upper_right' )
#  sz_dens.annotate_grids()
   sz_dens.set_buff_size( 1024 )
   sz_dens.save( mpl_kwargs={"dpi":300} )

   sz_real = yt.SlicePlot( ds, 'z', 'Real', center='c' )
#  sz_real.set_zlim( 'Real', 1.0e-1, 5.0e5 )
   sz_real.annotate_timestamp( corner='upper_right' )
   sz_real.set_buff_size( 1024 )
   sz_real.save( mpl_kwargs={"dpi":300} )

   sz_imag = yt.SlicePlot( ds, 'z', 'Imag', center='c' )
#  sz_imag.set_zlim( 'Imag', 1.0e-1, 5.0e5 )
   sz_imag.annotate_timestamp( corner='upper_right' )
   sz_imag.set_buff_size( 1024 )
   sz_imag.save( mpl_kwargs={"dpi":300} )
