# `configure.py` options
- Must enable
   - [[--model=ELBDM | Installation:-Option-List#--modle]]
   - [[--double | Installation:-Option-List#--double]]
   - [[--passive=1 | Installation:-Option-List#--passive]]
- Must disable
   - [[--gravity | Installation:-Option-List#--gravity]]
   - [[--particle | Installation:-Option-List#--particle]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Evolve the plane wave for six periods
2. Apply the periodic BC
   --> Set [[OPT__BC_FLU_* | Hydro#OPT__BC_FLU_XM]] = 1


# Note
1. Only support 1D --> Use `PWave_XYZ` to control the propagation direction
