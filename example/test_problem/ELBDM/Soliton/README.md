Compilation flags:
========================================
Enable : MODEL=ELBDM, GRAVITY
Disable: PARTICLE, UNSPLIT_GRAVITY, COMOVING


Default setup:
========================================
1. NX0_TOT_X/Y/Z = 96
2. MAX_LEVEL     = 2


Note:
========================================
1. Soliton_N = 1: test a single soliton
                  --> Set Soliton_RSeed < 0
                      Set Soliton_EmptyRegion = 0.0
   Soliton_N > 1: soliton merger
                  --> Set Soliton_RSeed >= 0
                      Set Soliton_EmptyRegion > 0.0 (e.g., 1.0e2)
