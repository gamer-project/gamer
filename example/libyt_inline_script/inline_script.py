import yt_libyt
import yt

yt.enable_parallelism()

def yt_inline():
    ds = yt_libyt.libytDataset()
    sz = yt.ProjectionPlot(ds, 'z', ('gamer', 'Temp'), center='c')

    if yt.is_root():
        sz.save()


def yt_inline_inputArg( fields ):
    ds = yt_libyt.libytDataset()
    sz = yt.ProjectionPlot(ds, 'z', fields, center='c')

    if yt.is_root():
        sz.save()
