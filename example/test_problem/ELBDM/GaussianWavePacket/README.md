# `configure.py` options
- Must enable
   - [[--model=ELBDM | Installation:-Option-List#--model]]
- Must disable
   - [[--gravity | Installation:-Option-List#--gravity]]
   - [[--particle | Installation:-Option-List#--particle]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Evolve the Gaussian wave packet for half of box
2. Apply the analytical solution as user-defined BC
   --> Set [[OPT__BC_FLU_* | Hydro#OPT__BC_FLU_XM]] = 4


# Note
1. Only support 1D --> Use `Gau_XYZ` to control the propagation direction
