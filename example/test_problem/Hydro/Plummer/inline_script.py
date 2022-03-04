import yt

yt.enable_parallelism()

def yt_inline():
    # Get data
    ds = yt.frontends.libyt.libytDataset()

    # Do ProjectionPlot to field Cloud0.
    sz = yt.ProjectionPlot(ds, 'z', ('gamer', 'Cloud0'), center='c')

    # Do ParticlePlot
    par = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y', 'particle_mass', center='c')

    if yt.is_root():
        sz.save()
        par.save()

def yt_inline_inputArg( fields ):
    # Get data
    ds = yt.frontends.libyt.libytDataset()

    # Do ProjectionPlot to fields input by gamer
    sz = yt.ProjectionPlot(ds, 'z', fields, center='c')

    if yt.is_root():
        sz.save()
