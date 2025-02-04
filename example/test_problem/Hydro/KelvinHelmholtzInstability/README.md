# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
- Must disable
  - [[--comoving | Installation:-Option-List#--comoving]]
  - [[--particle | Installation:-Option-List#--particle]]
  - [[--gravity | Installation:-Option-List#--gravity]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Refinement criteria: vorticity and regions near the shear plane
- Maximum refinement level ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]])=`3`


# Note
- Adopt the periodic boundary condition
- Shear velocity is defined on the x-y plane
- Random values are added to velocity in all three directions
- yt script `plot_density.py` for visualization
- Must disable `OPT__INIT_GRID_WITH_OMP`
  - Otherwise all threads would share the same random seed
