Compilation flags:
========================================
Enable :
        MODEL=HYDRO
        EOS=EOS_GAMMA_CR
        FLOAT8
        COSMIC_RAY
Disable:
        COMOVING
        PARTICLE
        GRAVITY
        CR_DIFFUSION


Default setup:
========================================
None


Note:
========================================
1. Support both 1D and 3D cases --> set by CR_Acoustic_Dir
2. A simple python script "plot_wave.py" is attached
3. "Record__L1Err" records the L1 errors
4. CosmicRay_AcousticWave_*_* record the numerical and analytical solutions along the diagonal
