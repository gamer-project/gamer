import yt

def yt_inline():
    ds = yt.frontends.libyt.libytDataset()
    sz = yt.SlicePlot( ds, 'z', 'Dens', center='c' )

#   sz.set_unit( 'Dens', 'msun/kpc**3' )
#   sz.set_zlim( 'Dens', 1.0e0, 1.0e6 )
    sz.annotate_grids( periodic=False )

    sz.save()

