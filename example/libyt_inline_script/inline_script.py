import yt

def yt_inline():
    ds = yt.frontends.libyt.libytDataset()
    sz = yt.SlicePlot( ds, 'z', 'Dens', center='c' )
    
    sz.set_zlim('Dens', 0.0, 3.5)
    sz.set_log('Dens', False)
    sz.set_cmap('Dens', 'arbre')
    sz.set_unit('Dens', 'code_mass/code_length**3')
    sz.set_axes_unit('code_length')
    sz.annotate_timestamp(time_unit='code_time', corner='upper_right', time_format='t = {time:.4f} {units}')
    sz.annotate_grids(periodic=False)
    sz.save(mpl_kwargs={"dpi":150})

    sz.save()

