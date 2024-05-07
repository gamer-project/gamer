import numpy as np
import yt

N = 11
ana_obj = 'CRay'
for i in range(N):
    ds = yt.load('./Data_%06d'%i)
    p = yt.SlicePlot(ds, 'y', ana_obj, center=[0.5, 0.5, 0.5])
    p.set_log(ana_obj, False)
    p.set_zlim(ana_obj, zmin=0.0, zmax=1.5)
    p.annotate_grids()
    p.save()
