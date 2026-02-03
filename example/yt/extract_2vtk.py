"""
Purpose:
    Extracts uniform 3D data from yt datasets and exports them to VTK format

Usage:
    1. Install the requriments via `pip install yt h5py pyevtk`
    2. Set up the "fn_in" and "fn_out"
    3. Adjust the "left_edge" and "right_edge" to the target cube boundaries
    4. Adjust the "dims" to specify the number of cells in each dimension
    5. Execute the script via `python extract_2vtk.py`

Requirements:
    - yt
    - pyevtk
    - h5py
"""

import yt
from pyevtk.hl import imageToVTK


# (1) Set up the input and output files
#     Note: no extension needed in "fn_out"
fn_in  = "Data_000000"
fn_out = "Data_000000_VTK"


# (2) Extract data using yt
ds = yt.load(fn_in)

# (2.1) Define the spatial extent and resolution of the uniform grid
center     = ds.domain_center
width      = ds.quan(100, "code_length")
left_edge  = center - 0.5 * width
right_edge = center + 0.5 * width
dims       = [128, 128, 128]

# (2.2) Sample data onto a uniform grid using arbitrary_grid()
ag = ds.arbitrary_grid(left_edge = left_edge, right_edge = right_edge, dims = dims)

# (2.3) Extract fields and convert to CGS units
dens  = ag["gas", "density"   ].in_cgs().v
vel_x = ag["gas", "velocity_x"].in_cgs().v
vel_y = ag["gas", "velocity_y"].in_cgs().v
vel_z = ag["gas", "velocity_z"].in_cgs().v


# (3) Export cell-centered data to VTK
#     Note: The origin and spacing are in code units
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
