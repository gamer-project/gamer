# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
   - [[MHD | Installation: Simulation-Options#MHD]]
- Must disable
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. Periodic BC
2. Refinement criterion: current density ([[OPT__FLAG_CURRENT | Runtime-Parameters:-Refinement#OPT__FLAG_CURRENT]])
3. [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]] = 2


# Note
1. Reference: [Orszag & Tang, 1998, J. Fluid Mech., 90, 129](https://doi.org/10.1017/S002211207900210X)
2. No z dependence
