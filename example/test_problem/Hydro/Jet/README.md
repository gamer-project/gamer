# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]] (when enabling `JET_HSE`)
- Must disable
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]] (when disabling `JET_HSE`)
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
None


# Note
1. Ref: [Molnar, S., Schive, H-Y., et al., 2017, ApJ, 835, 57](https://arxiv.org/abs/1612.02341)
2. Recommended boundary conditions
   - With `JET_HSE` and large `Jet_BgVel_y`: user-defined BC on -y and outflow BC on other faces
   - Without `JET_HSE`: outflow BC on all faces
