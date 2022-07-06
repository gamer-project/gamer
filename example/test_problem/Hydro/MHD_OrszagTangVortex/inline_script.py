import yt

yt.enable_parallelism()

def yt_inline():
    ds = yt.frontends.libyt.libytDataset()
    slc = yt.OffAxisSlicePlot(ds, [1, 1, 0], [("gas", "density")], center="c")
    slc.annotate_cquiver(("gas", "cutting_plane_velocity_x"), ("gas", "cutting_plane_velocity_y"), factor=10, plot_args={"color":"orange"}, )

    # Currently, saving annotate_cquiver will only work when it is outside of yt.is_root() clause like this.
    slc.save()

def yt_inline_inputArg( field ):
    pass

