"""
This file generates a toy vector potential for import into GAMER using the
OPT__INIT_BFIELD_BYFILE parameter. It does the following:

1. Generates a uniform coordinate grid
2. Defines a vector potential on the coordinate grid
3. Saves the coordinate grid and the vector potential to an HDF5 file

The units of the vector potential and the coordinate arrays should be the same
as those used in GAMER. So:

* coordinates are in UNIT_L
* vector potential components are in UNIT_B*UNIT_L = sqrt(4*pi*UNIT_P)*UNIT_L
  where UNIT_P = UNIT_M/UNIT_L/UNIT_T**2

The file should also be named "B_IC" for GAMER to recognize it.

It requires NumPy, h5py, and HDF5 to be installed.
"""

import h5py
import numpy as np

# Number of cells along each dimension of the input grid.
# This is somewhat arbitrary, but should be chosen in
# such a way as to adequately resolve the vector potential.

ddims = np.array([128]*3, dtype='int')

# Left edge and right edge coordinates of the desired
# simulation domain which will be used in GAMER.

le = np.zeros(3)
re = np.ones(3)

# Since we need to take derivatives of the vector potential
# to get the magnetic field on the simulation domain, the
# input grid must be extended a bit beyond this boundary.
# We therefore add a buffer of three cells on each side.
# (Three cells are necessary to solve some corner cases
# resulting from round-off errors.)

delta = (re-le)/ddims
ddims += 6
le -= 3.0*delta
re += 3.0*delta

# Construct the grid cell edge coordinates

x = np.linspace(le[0], re[0], ddims[0]+1)
y = np.linspace(le[1], re[1], ddims[1]+1)
z = np.linspace(le[2], re[2], ddims[2]+1)

# Find the grid cell midpoints

x = 0.5*(x[1:]+x[:-1])
y = 0.5*(y[1:]+y[:-1])
z = 0.5*(z[1:]+z[:-1])

# Use the 1-D coordinate arrays to consruct 3D coordinate arrays
# that we will use to compute an analytic vector potential

xx, yy, zz = np.meshgrid(x, y, z, sparse=False, indexing='ij')

# Toy vector potential which depends on all three coordinates

Ax = 3.0*yy*zz*zz
Ay = 2.0*xx*xx*zz
Az = yy*yy*xx

# Write the ICs to an HDF5 file

f = h5py.File("B_IC", "w")

# Write coordinate arrays

f.create_dataset("x", data=x)
f.create_dataset("y", data=y)
f.create_dataset("z", data=z)

#  Write vector potential arrays

f.create_dataset("magnetic_vector_potential_x", data=Ax)
f.create_dataset("magnetic_vector_potential_y", data=Ay)
f.create_dataset("magnetic_vector_potential_z", data=Az)

# Close the file

f.flush()
f.close()
