Compilation flags:
========================================
Enable : MODEL=ELBDM, GRAVITY, COMOVING
Disable: PARTICLE


Default setup:
========================================
1. ELBDM_MASS              8.0e-23 (eV/c^2)
   A_INIT                  3.124024e-4
   OMEGA_M0                0.2835
   HUBBLE0                 0.6955
   BOX_SIZE                1.4 (Mpc/h)

2. MAX_LEVEL               3
   OPT__FLAG_RHO           1
   OPT__FLAG_LOHNER_DENS   1
   --> Input__Flag_Rho    = example/input/Input__Flag_Rho8
       Input__Flag_Lohner = example/input/Input__Flag_Lohner__FLASH1_LSS_DH3.0


Note:
========================================
1. Cosmological large-scale structure simulations

2. Download the IC file: sh download_ic.sh

3. Default maximum resolution is low
   --> Only 1.4/(256*2^3) ~ 0.7 kpc/h (comoving)
   --> Can only marginally resolve a central soliton

4. Some yt visualization scripts are put in "plot_script"
