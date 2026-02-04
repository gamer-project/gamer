# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
  - [[--gravity | Installation:-Option-List#--gravity]] (when enabling `JET_HSE`)
- Must disable
  - [[--gravity | Installation:-Option-List#--gravity]] (when disabling `JET_HSE`)
  - [[--particle | Installation:-Option-List#--particle]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
None


# Note
- Ref: [Molnar, S., Schive, H-Y., et al., 2017, ApJ, 835, 57](https://arxiv.org/abs/1612.02341)
- Recommended boundary conditions
  - With `JET_HSE` and large `Jet_BgVel_y`: user-defined BC on -y and outflow BC on other faces
  - Without `JET_HSE`: outflow BC on all faces
