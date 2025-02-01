Compilation flags:
========================================
Enable : MODEL=HYDRO, GRAVITY (when enabling JET_HSE)
Disable: GRAVITY (when disabling JET_HSE), PARTICLE


Default setup:
========================================


Note:
========================================
1. Ref: Molnar, S., Schive, H-Y., et al., 2017, ApJ, 835, 57
        https://arxiv.org/abs/1612.02341
2. Recommended boundary conditions:
      With JET_HSE and large Jet_BgVel_y: user-defined BC on -y and outflow BC on other faces
      Without JET_HSET: outflow BC on all faces
