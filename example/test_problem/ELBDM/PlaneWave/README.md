# Compilation flags
- Must enable
   - [[MODEL=ELBDM | Installation: Simulation-Options#MODEL]]
   - [[FLOAT8 | Installation: Simulation-Options#FLOAT8]]
   - [[NCOMP_PASSIVE_USER=1 | Installation: Simulation-Options#NCOMP_PASSIVE_USER]]
- Must disable
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. Evolve the plane wave for six periods
2. Apply the periodic BC
   --> Set [[OPT__BC_FLU_* | Hydro#OPT__BC_FLU_XM]] = 1


# Note
1. Only support 1D --> Use `PWave_XYZ` to control the propagation direction
