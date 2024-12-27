import yt_libyt
import yt
import numpy as np

yt.enable_parallelism()

def yt_inline():
# ref: https://yt-project.org/doc/examining/low_level_inspection.html#examining-grid-data-in-a-fixed-resolution-array
   ds = yt_libyt.libytDataset()
   idx = int(str(ds).split("_")[0].split("g")[-1])

   # target field
   field = "Dens"

   # target amr level
   lv = 1

   # get unsmoothed/smoothed fixed resolution array
   ad = ds.covering_grid( level=lv, left_edge=[0.04375,0.04375,0.04375], dims=512 )

   density = ad[field]

   np.savez("covering-grid_test_Data_%06d_lv=%d.npz"%(idx,lv), Dens=np.array(density).astype(np.float32))

def yt_inline_inputArg( fields ):
   pass
