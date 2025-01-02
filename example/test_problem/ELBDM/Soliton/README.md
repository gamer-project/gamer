# `configure.py` options
- Must enable
   - [[--model=ELBDM | Installation:-Option-List#--modle]]
   - [[--gravity | Installation:-Option-List#--gravity]]
- Must disable
   - [[--comoving | Installation:-Option-List#--comoving]]
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--unsplit_gravity | Installation:-Option-List#--unsplit_gravity]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. [[NX0_TOT_X/Y/Z | Runtime-Parameters:-General#NX0_TOT_X]] = 96
2. [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]] = 2


# Note
1. `Soliton_N` = 1: test a single soliton
   * Set `Soliton_RSeed` < 0
   * Set `Soliton_EmptyRegion` = 0.0

2. `Soliton_N` > 1: soliton merger
   * Set `Soliton_RSeed` >= 0
   * Set `Soliton_EmptyRegion` > 0.0 (e.g., 1.0e2)
