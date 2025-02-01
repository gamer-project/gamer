Compilation flags:
========================================
Enable : MODEL=HYDRO
Disable: COMOVING, PARTICLE, GRAVITY


Default setup:
========================================
1. Adopt the Lohner's error estimator of both mass density and pressure as the refinement criteria
2. Refinement threshold of the Lohner's error estimator = 0.80
3. Maximum refinement level (MAX_LEVEL) = 2


Note:
========================================
1. Ref: Ryu, D., Ostriker, J. P., Kang, H., & Cen, R. 1993, ApJ, 414, 1

2. This test is good for testing the dual-energy formalism
   --> Without it the preshock region will be over-heated

