#!/usr/bin/env python

import yt
import glob

fns = glob.glob("Data_000*[0-9]")
fns.sort()

ts = yt.DatasetSeries(fns)

for ds in ts:
    slc = yt.SlicePlot(ds, "z", ["velocity_magnitude","particle_density_on_grid"])
    slc.annotate_particles(1.0, p_size=10)
    slc.annotate_velocity(plot_args={"color": "red"})
    slc.annotate_grids()
    slc.save()

