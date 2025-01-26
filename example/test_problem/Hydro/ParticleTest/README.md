# `configure.py` options
- Must enable
   - [[--model | Installation:-Option-List#--model]]=`HYDRO`
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--gravity | Installation:-Option-List#--gravity]]
   - [[--tracer | Installation:-Option-List#--tracer]]
   - [[--hdf5 | Installation:-Option-List#--hdf5]]
- Must disable
   - [[--comoving | Installation:-Option-List#--comoving]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Maximum refinement level ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]])=`3`
2. Both massive and tracer particles enabled

# Note
1. Test particle capabilities: 2-body circular orbit for massive particles, a
   lattice of tracer particles which rotates due to a circular velocity field
   from the fluid with constant angular frequency.
2. Set `ParTest_Use_Tracers = 1` to use tracer particles, `ParTest_Use_Massive = 1`
   to use massive particles. At least one must be enabled.
3. No need to set [[PAR_NPAR | Runtime-Parameters:-Particles#PAR_NPAR]] in `Input__Parameter`
   since it will be reset automatically.
4. [[OPT__FREEZE_FLUID | Runtime-Parameters:-Hydro#OPT__FREEZE_FLUID]] must be 1 for this test.
5. If [[OPT__FLAG_USER | Runtime-Parameters:-Refinement#OPT__FLAG_USER]] is set,
   a rectangle will be refined in the center of the domain.
> [!NOTE]
> The active particle orbit is _not_ stable in this case due to self-force at refinement boundaries.
6. Density of the gas must be kept very small to avoid gravitational effects from the gas.
7. Use python script `make_tracer_pngs.py` with yt to make PNG images, use
   bash script `make_tracer_movie.sh` to combine PNGs into movies with `ffmpeg`.
8. Set [[OPT__OUTPUT_PAR_MESH | Runtime-Parameters:-Outputs#OPT__OUTPUT_PAR_MESH]]=`1` to output the mesh attributes, as listed in
   `Input__Par_Mesh`, to tracer particles. Use the python script `check_mesh2tracer.py`
   to verify the mesh attributes recorded in the HDF5 snapshots.