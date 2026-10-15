import yt_libyt
import yt

yt.enable_parallelism()

def yt_inline():
    # Get data
    ds = yt_libyt.libytDataset()

    # Do ProjectionPlot to density field
    sz = yt.ProjectionPlot(ds, 'z', ('gamer', 'Dens'), center='c')

    # Do ParticlePlot
    par = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y', 'particle_mass', center='c')

    # Do ProjectionPlot to field Cloud0.
    if ('gamer', 'Cloud0') in ds.field_list:
        scz = yt.ProjectionPlot(ds, 'z', ('gamer', 'Cloud0'), center='c')

    if yt.is_root():
        sz.save()
        par.save()
        if ('gamer', 'Cloud0') in ds.field_list:
            scz.save()

def yt_inline_inputArg( fields ):
    pass
