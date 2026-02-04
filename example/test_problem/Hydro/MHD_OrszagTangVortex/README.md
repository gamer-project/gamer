# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
  - [[--mhd | Installation:-Option-List#--mhd]]
- Must disable
  - [[--particle | Installation:-Option-List#--particle]]
  - [[--gravity | Installation:-Option-List#--gravity]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Periodic BC
- Refinement criterion: current density ([[OPT__FLAG_CURRENT | Runtime-Parameters:-Refinement#OPT__FLAG_CURRENT]])
- [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]=`2`


# Note
- Reference
  - [Orszag & Tang, 1998, J. Fluid Mech., 90, 129](https://doi.org/10.1017/S002211207900210X)
- No z dependence
