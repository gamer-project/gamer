# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
   - [[FLOAT8 | Installation: Simulation-Options#FLOAT8]]
- Must disable
   - [[COMOVING | Installation: Simulation-Options#COMOVING]]
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


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
      - [[SRHD | Installation: Simulation-Options#SRHD]]
      - [[EOS=EOS_TAUBMATHEWS | Installation: Simulation-Options#EOS]]
      - [[FLU_SCHEME=MHM | Installation: Simulation-Options#FLU_SCHEME]]
      - [[RSOLVER=HLLC | Installation: Simulation-Options#RSOLVER]]
      - [[LR_SCHEME=PLM | Installation: Simulation-Options#LR_SCHEME]]
   - In `Input__TestProb`,
      - `Acoustic_v0` and `Acoustic_Cs` will be useless
      - set `Acoustic_Temp_Bg = 1.0e+10` for high-temperature case
      - set `Acoustic_Temp_Bg = 1.0e-10` for low-temperature case
   - To check the L1 error convergence rate
      - A L1 error plotting script `plot_L1error_SRHD.py` is attached
      - An automatic script `Run_L1ErrorConvergenceTest_SRHD.sh` is attached
      - Run `sh Run_L1ErrorConvergenceTest_SRHD.sh`
      - Reference: Figure 3. in [Tseng et al. 2021, MNRAS, 504, 3298](https://doi.org/10.1093/mnras/stab1006)