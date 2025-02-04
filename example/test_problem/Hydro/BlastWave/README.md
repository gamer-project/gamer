# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
- Must disable
  - [[--particle | Installation:-Option-List#--particle]]
  - [[--gravity | Installation:-Option-List#--gravity]]
- Available options
  - [[--mhd | Installation:-Option-List#--mhd]]
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Adopt the Lohner's error estimator on pressure as the refinement criteria
  - `Refinement threshold`=`0.80`

- Maximum refinement level ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]])=`3`


# Note
- Two plot scripts are attached: `plot_profile.gpt` and `plot_slice.py`
- To reset the magnetic field during the runtime, turn on the following options.
  - [[OPT__RESET_FLUID | Runtime-Parameters:-Hydro#OPT__RESET_FLUID]]

  - `Blast_ResetB_VecPot` (optional but recommended)

    Improve divergence-free

  - [[OPT__FLAG_USER | Runtime-Parameters:-Refinement#OPT__FLAG_USER]] (optional but recommended)

    Refine regions within `~ 6*Blast_ResetB_r0` to the maximum level to improve divergence-free
    (since even using vector potential, divergence-free can still break at the coarse-fine boundaries)
