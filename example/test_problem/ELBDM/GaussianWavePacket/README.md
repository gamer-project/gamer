# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`ELBDM`
- Must disable
  - [[--gravity | Installation:-Option-List#--gravity]]
  - [[--particle | Installation:-Option-List#--particle]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Evolve the Gaussian wave packet for half of box
- Apply the analytical solution as user-defined BC
  - Set [[OPT__BC_FLU_* | Runtime-Parameters:-Hydro#OPT__BC_FLU_XM]]=`4`


# Note
- Only support 1D
  - Use `Gau_XYZ` to control the propagation direction
