# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
  - [[--comoving | Installation:-Option-List#--comoving]]
  - [[--gravity | Installation:-Option-List#--gravity]]
- Must disable
  - [[--particle | Installation:-Option-List#--particle]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Refinement criteria: density magnitude and Lohner's error estimator for density
  - Refinement threshold of the Lohner error estimator = 0.80
- [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]=`2`


# Note
- Dual-energy formalism is important for evolving pressure accurately in the preshock region
