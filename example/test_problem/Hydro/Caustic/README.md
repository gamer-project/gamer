# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
- Must disable
   - [[COMOVING | Installation: Simulation-Options#COMOVING]]
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. Adopt the Lohner's error estimator of both mass density and pressure as the refinement criteria
2. Refinement threshold of the Lohner's error estimator = 0.80
3. Maximum refinement level ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]) = 2


# Note
1. Ref: [Ryu, D., Ostriker, J. P., Kang, H., & Cen, R. 1993, ApJ, 414, 1](https://doi.org/10.1086/173051)
2. This test is good for testing the dual-energy formalism
   --> Without it the preshock region will be over-heated

