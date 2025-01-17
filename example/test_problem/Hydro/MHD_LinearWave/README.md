# `configure.py` options
- Must enable
   - [[--model | Installation:-Option-List#--model]]
   - [[--mhd | Installation:-Option-List#--mhd]]
   - [[--double | Installation:-Option-List#--double]]
- Must disable
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--gravity | Installation:-Option-List#--gravity]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Resolution = 64^3
2. Run for one period


# Note
1. Support both 1D and 3D cases --> set by `MHDLinear_Dir`
2. A simple gnuplot script `plot.gpt` is attached
3. `Record__L1Err` records the L1 errors
4. `MHDLinearWave_*_*` record the numerical and analytical solutions along the diagonal
