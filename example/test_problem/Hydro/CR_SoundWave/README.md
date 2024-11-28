# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
   - [[EOS=EOS_GAMMA_CR | Installation: Simulation-Options#EOS]]
   - [[FLOAT8 | Installation: Simulation-Options#FLOAT8]]
   - [[COSMIC_RAY | Installation: Simulation-Options#COSMIC_RAY]]
- Must disable
   - [[COMOVING | Installation: Simulation-Options#COMOVING]]
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
   - [[CR_DIFFUSION | Installation: Simulation-Options#CR_DIFFUSION]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. Resolution = 64^3
2. Run for one period


# Note
1. Support both 1D and 3D cases --> set by `CR_Acoustic_Dir`
2. A simple Python script `plot_wave.py` is attached
3. `Record__L1Err` records the L1 errors
4. `CosmicRay_AcousticWave_*_*` record the numerical and analytical solutions along the diagonal
