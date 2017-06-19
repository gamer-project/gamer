import yt

yt.enable_parallelism()

idx_start = 0
idx_end   = 8
didx      = 1
prefix    = '../Data_'

ts = yt.load( [ prefix+'%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

#  density
   sz_dens = yt.SlicePlot( ds, 'z', 'density', center='c' )

   sz_dens.set_zlim( 'density', 1.0e-30, 1.0e-25 )
   sz_dens.set_cmap('density', 'algae')
#  sz_dens.set_figure_size(16)
#  sz_dens.set_buff_size(2048)
#  sz_dens.set_font( {'size':24} )
   sz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_dens.annotate_grids( periodic=False )

   sz_dens.save()


#  temperature
   sz_temp = yt.SlicePlot( ds, 'z', 'temperature', center='c' )

   sz_temp.set_zlim( 'temperature', 2.0e7, 1.5e8 )
   sz_temp.set_cmap('temperature', 'afmhot')
#  sz_temp.set_figure_size(16)
#  sz_temp.set_buff_size(2048)
#  sz_temp.set_font( {'size':24} )
   sz_temp.annotate_timestamp( time_unit='Myr', corner='upper_right' )
#  sz_temp.annotate_grids( periodic=False )

   sz_temp.save()
