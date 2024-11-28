# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
- Must disable
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
- Available options
   - [[MHD | Installation: Simulation-Options#MHD]]
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. Adopt the Lohner's error estimator on pressure as the refinement criteria
   --> `Refinement threshold = 0.80`
2. Maximum refinement level ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]] = 3)


# Note
1. Two plot scripts are attached: `plot_profile.gpt` and `plot_slice.py`
2. To reset the magnetic field during the runtime, turn on the following options.
   1. [[OPT__RESET_FLUID | Hydro#OPT__RESET_FLUID ]]
   2.  `Blast_ResetB_VecPot` (optional but recommended)

       Improve divergence-free
   3.  [[OPT__FLAG_USER | Runtime-Parameters:-Refinement#OPT__FLAG_USER ]] (optional but recommended)

       Refine regions within `~ 6*Blast_ResetB_r0` to the maximum level to improve divergence-free
       (since even using vector potential, divergence-free can still break at the coarse-fine boundaries)
