# `configure.py` options
- Must enable
   - [[--model | Installation:-Option-List#--modle]]
   - [[--double | Installation:-Option-List#--double]]
- Must disable
   - [[--comoving | Installation:-Option-List#--comoving]]
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--gravity | Installation:-Option-List#--gravity]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Resolution = 64^3
2. Run for one period


# Note
1. Support both 1D and 3D cases --> set by `Acoustic_Dir`
2. A simple gnuplot script `plot.gpt` is attached
3. `Record__L1Err` records the L1 errors
4. `AcousticWave_*_*` record the numerical and analytical solutions along the diagonal
5. For SRHD:
   - Also enable these compilation flags:
      - [[--srhd | Installation:-Option-List#--srhd]]
      - [[--eos=TAUBMATHEWS | Installation:-Option-List#--eos]]
      - [[--flu_scheme=MHM | Installation:-Option-List#--flu_scheme]]
      - [[--flux=HLLC | Installation:-Option-List#--flux]]
      - [[--slope=PLM | Installation:-Option-List#--slope]]
   - In `Input__TestProb`,
      - `Acoustic_v0` and `Acoustic_Cs` will be useless
      - set `Acoustic_Temp_Bg = 1.0e+10` for high-temperature case
      - set `Acoustic_Temp_Bg = 1.0e-10` for low-temperature case
   - To check the L1 error convergence rate
      - A L1 error plotting script `plot_L1error_SRHD.py` is attached
      - An automatic script `Run_L1ErrorConvergenceTest_SRHD.sh` is attached
      - Run `sh Run_L1ErrorConvergenceTest_SRHD.sh`
      - Reference: Figure 3. in [Tseng et al. 2021, MNRAS, 504, 3298](https://doi.org/10.1093/mnras/stab1006)
