# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
- Must disable
  - [[--comoving | Installation:-Option-List#--comoving]]
  - [[--particle | Installation:-Option-List#--particle]]
  - [[--gravity | Installation:-Option-List#--gravity]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Adopt the Lohner's error estimator of both mass density and pressure as the refinement criteria
- Refinement threshold of the Lohner's error estimator = 0.80
- Maximum refinement level ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]])=`2`


# Note
- Ref: [Ryu, D., Ostriker, J. P., Kang, H., & Cen, R. 1993, ApJ, 414, 1](https://doi.org/10.1086/173051)
- This test is good for testing the dual-energy formalism
  - Without it the preshock region will be over-heated
