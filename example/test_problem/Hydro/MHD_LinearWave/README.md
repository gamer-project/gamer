# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
   - [[MHD | Installation: Simulation-Options#MHD]]
   - [[FLOAT8 | Installation: Simulation-Options#FLOAT8]]
- Must disable
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. Resolution = 64^3
2. Run for one period


# Note
1. Support both 1D and 3D cases --> set by `MHDLinear_Dir`
2. A simple gnuplot script `plot.gpt` is attached
3. `Record__L1Err` records the L1 errors
4. `MHDLinearWave_*_*` record the numerical and analytical solutions along the diagonal
