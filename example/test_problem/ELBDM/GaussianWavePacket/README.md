# Compilation flags
- Must enable
   - [[MODEL=ELBDM | Installation: Simulation-Options#MODEL]]
- Must disable
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. Evolve the Gaussian wave packet for half of box
2. Apply the analytical solution as user-defined BC
   --> Set [[OPT__BC_FLU_* | Hydro#OPT__BC_FLU_XM]] = 4


# Note
1. Only support 1D --> Use `Gau_XYZ` to control the propagation direction
