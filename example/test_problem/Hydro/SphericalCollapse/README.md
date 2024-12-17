# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
   - [[COMOVING | Installation: Simulation-Options#COMOVING]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
- Must disable
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. Refinement criteria: density magnitude and Lohner's error estimator for density
   --> Refinement threshold of the Lohner error estimator = 0.80
2. [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]] = 2


# Note
1. Dual-energy formalism is important for evolving pressure accurately in the
   preshock region

