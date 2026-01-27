# Extract uniform 3D data using yt and export to VTK files

import yt
from pyevtk.hl import imageToVTK


# setting
# --> no extension needed in "fn_out"
fn_in  = "Data_000000"
fn_out = "Data_000000_VTK"


# use arbitrary_grid() to sample the data on a uniform grid
ds = yt.load(fn_in)

center     = ds.domain_center
width      = ds.quan(100, "code_length")
left_edge  = center - 0.5 * width
right_edge = center + 0.5 * width
dims       = [128, 128, 128]

ag = ds.arbitrary_grid(left_edge = left_edge, right_edge = right_edge, dims = dims)

dens  = ag["gas", "density"   ].in_cgs().v
vel_x = ag["gas", "velocity_x"].in_cgs().v
vel_y = ag["gas", "velocity_y"].in_cgs().v
vel_z = ag["gas", "velocity_z"].in_cgs().v


# export cell centered data to VTK
origin  = tuple(left_edge.in_units("code_length").v)
spacing = tuple(width.v / dim  for dim in dims)

imageToVTK(
    fn_out,
    origin   = origin,
    spacing  = spacing,
    cellData = {
        "density" : dens,
        "velocity": (vel_x, vel_y, vel_z)
    }
)
