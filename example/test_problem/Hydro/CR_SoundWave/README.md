# `configure.py` options
- Must enable
   - [[--model | Installation:-Option-List#--model]]
   - [[--eos=COSMIC_RAY | Installation:-Option-List#--eos]]
   - [[--double | Installation:-Option-List#--double]]
   - [[--cosmic_ray | Installation:-Option-List#--cosmic_ray]]
- Must disable
   - [[--comoving | Installation:-Option-List#--comoving]]
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--gravity | Installation:-Option-List#--gravity]]
   - [[--cr_diffusion | Installation:-Option-List#--cr_diffusion]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Resolution = 64^3
2. Run for one period


# Note
1. Support both 1D and 3D cases --> set by `CR_Acoustic_Dir`
2. A simple Python script `plot_wave.py` is attached
3. `Record__L1Err` records the L1 errors
4. `CosmicRay_AcousticWave_*_*` record the numerical and analytical solutions along the diagonal
