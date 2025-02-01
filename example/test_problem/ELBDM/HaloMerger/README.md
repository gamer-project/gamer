Compilation flags:
========================================
Enable : MODEL=ELBDM, GRAVITY (, PARTICLE, SUPPORT_GSL)
Disable: COMOVING


Default setup:
========================================
1. Code units
   (1) UNIT_L = Mpc/h, where h=0.6955 is the present dimensionless Hubble parameter
              = 1437.814521231748 kpc
   (2) UNIT_V = 100 km/s
   (3) UNIT_D = rho_bg (background matter density at z=0)
              = 38.06 Msun/kpc^3
       --> Mass density and wavefunction are normalized to rho_bg
   (4) UNIT_T = 14068.4678922741 Myr
   (5) UNIT_M = 1.131e+11 Msun

2. ELBDM_MASS              1.0e-22
   ELBDM_REMOVE_MOTION_CM  0
   ELBDM_TAYLOR3_AUTO      0

3. OPT__BC_FLU_*           1  (periodic)
   OPT__BC_POT             2  (isolated)

4. MAX_LEVEL               3
   OPT__FLAG_RHO           1

5. END_T                   0.25

Note:
========================================
1. Simulate the merger of halos and solitons in an external potential

2. Edit Input_TestProb_Halo to specify the parameters for the halos
   a. Add new parameters with increasing indexes if the number of halos is more than 2
   b. For HaloMerger_Halo_InitMode == 1
      - The initial condition of halos is constructed by reading the HALO_IC file as a table and performing linear interpolation
         - Note that this HALO_IC has nothing to do with the built-in option OPT__INIT = 3
      - The HALO_IC must be single AMR level
         - If there is a halo UM_IC (for OPT__INIT = 3) that has multiple AMR levels (and there is Input__UM_IC_RefineRegion),
           the Python script Make_UM_IC_uniform.py is provided to convert it to
           a single-level UM_IC with the specified level, which can be used as the HALO_IC
           (However, those levels higher than Target_lv+log_2( 2*PatchSize ) cannot be handled and will be ignored)
         - If converting the halo UM_IC to be single-level is not feasible,
           switching to the initialization option OPT__INIT = 3 to load a multi-level UM_IC is also an alternative
           (However, only one UM_IC can be loaded and all the parameters in Input_TestProb_Halo will not be used)
      - Note that HaloMerger_Halo_*_CenCoord* in Input__TestProb_Halo is the HALO_IC box center rather than the exact halo center

3. Edit Input_TestProb_Soliton to specify the parameters for the solitons
   a. Add new parameters with increasing indexes if the number of solitons is more than 2
   b. For HaloMerger_Soliton_InitMode == 1
      - The initial condition of solitons is constructed by reading the table of soliton density profile
        (the density profile will be rescaled to the given core radius or core density if HaloMerger_Soliton_*_DensProf_Rescale == 1)
   c. For HaloMerger_Soliton_InitMode == 2
      - The initial condition of solitons is constructed by using the analytical formula of soliton density profile

4. Edit Input_TestProb_ParCloud to specify the parameters for the particle clouds
   a. Add new parameters with increasing indexes if the number of particle clouds is more than 2
   b. For HaloMerger_ParCloud_InitMode == 1
      - The initial condition of particle clouds is constructed by reading the table of density profile and using Par_EquilibriumIC()
      - The default particle clouds use the HaloDensityProfile (see Note 7. below) to represent the CDM halos
      - Enable compilation options: PARTICLE, SUPPORT_GSL
      - Set the parameter OPT__FREEZE_FLUID to 1 in Input__Parameter for particle-only simulations

5. Turn on OPT__EXT_POT == 1 to add external potential
   a. The external potential of a uniform-density sphere, which is proportional to r^2 inside and proportional to 1/r outside

6. Download the default HALO_IC file: sh download_ic.sh
   a. Halo mass = 4.0960e+09 Msun
   b. Without soliton
   c. N = 640, L = 0.0646921095 Mpc/h, single-precision

7. Generate the default HaloDensityProfile and SolitonDensityProfile: python Make_DensityProfile.py
   a. The parameters for the halo density profile are used to fit the ELBDM halo

8. Some examples of yt visualization scripts are put in "plot_script"

9. The corresponding wavelength is 0.00083778702 Mpc/h when ELBDM_MASS = 1.0e-22 and velocity = 1.0*100 km/s
   -> Make sure the resolution is high enough to resolve the wavelength

10. The input halos and solitons may overlap with each other initially
   -> Their wavefunctions, instead of densities, are added directly
   -> Note that the interference will cause the density distribution to be different from the sum of individual density
