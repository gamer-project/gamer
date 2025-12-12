# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
  - [[--mhd | Installation:-Option-List#--mhd]]
  - [[--double | Installation:-Option-List#--double]]
- Must disable
  - [[--particle | Installation:-Option-List#--particle]]
  - [[--gravity | Installation:-Option-List#--gravity]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Resolution = 64^3
- Run for one period


# Note
- Support both 1D and 3D cases
  - set by `MHDLinear_Dir`
- A simple gnuplot script `plot.gpt` is attached
- `Record__L1Err` records the L1 errors
- `MHDLinearWave_*_*` record the numerical and analytical solutions along the diagonal
