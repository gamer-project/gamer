# Compilation flags
- Must enable
   - [[MODEL=ELBDM | Installation: Simulation-Options#MODEL]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
- Must disable
   - [[COMOVING | Installation: Simulation-Options#COMOVING]]
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[UNSPLIT_GRAVITY | Installation: Simulation-Options#UNSPLIT_GRAVITY]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


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
