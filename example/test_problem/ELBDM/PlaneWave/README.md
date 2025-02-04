# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`ELBDM`
  - [[--double | Installation:-Option-List#--double]]
  - [[--passive | Installation:-Option-List#--passive]]=`1`
- Must disable
  - [[--gravity | Installation:-Option-List#--gravity]]
  - [[--particle | Installation:-Option-List#--particle]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Evolve the plane wave for six periods
- Apply the periodic BC
  - Set [[OPT__BC_FLU_* | Runtime-Parameters:-Hydro#OPT__BC_FLU_XM]]=`1`


# Note
- Only support 1D
  - Use `PWave_XYZ` to control the propagation direction
