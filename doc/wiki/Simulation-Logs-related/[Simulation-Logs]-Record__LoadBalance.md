This file records the load-balancing information on each level in each
MPI process as well as the overall load imbalance. The load imbalance factor
is defined as <var>(Max-Ave)/Ave</var>, where <var>Max</var> and <var>Ave</var>
are the maximum workload of one process and the average workload among all processes,
respectively.

Example:
``` markdown
Time 6.8000000e-01,  Step     461,  NPatch       333632,  NPar      9923891

Rank           Level 0            Level 1            Level 2            Level 3            Level 4
   0 5.34e+03(  -0.33%) 7.52e+03(  -0.11%) 1.66e+04(  +0.28%) 2.29e+04(  +0.27%) 8.88e+04(  -0.22%)
   1 5.37e+03(  +0.30%) 7.53e+03(  +0.08%) 1.66e+04(  +0.33%) 2.28e+04(  -0.13%) 9.06e+04(  +1.83%)
   2 5.35e+03(  -0.16%) 7.51e+03(  -0.18%) 1.65e+04(  +0.08%) 2.28e+04(  -0.24%) 8.81e+04(  -0.96%)
   3 5.36e+03(  +0.07%) 7.54e+03(  +0.13%) 1.64e+04(  -0.57%) 2.28e+04(  -0.27%) 8.76e+04(  -1.54%)
   4 5.34e+03(  -0.27%) 7.53e+03(  +0.09%) 1.65e+04(  -0.01%) 2.28e+04(  -0.10%) 8.89e+04(  -0.15%)
   5 5.38e+03(  +0.44%) 7.53e+03(  +0.03%) 1.65e+04(  -0.22%) 2.28e+04(  +0.04%) 8.93e+04(  +0.33%)
   6 5.37e+03(  +0.20%) 7.52e+03(  -0.06%) 1.64e+04(  -0.68%) 2.27e+04(  -0.41%) 8.93e+04(  +0.38%)
   7 5.34e+03(  -0.26%) 7.53e+03(  +0.01%) 1.67e+04(  +0.80%) 2.30e+04(  +0.84%) 8.93e+04(  +0.33%)
---------------------------------------------------------------------------------------------------
Sum: 4.29e+04           6.02e+04           1.32e+05           1.83e+05           7.12e+05
Ave: 5.36e+03           7.53e+03           1.65e+04           2.28e+04           8.90e+04
Max: 5.38e+03           7.54e+03           1.67e+04           2.30e+04           9.06e+04
Imb:    0.44%              0.13%              0.80%              0.84%              1.83%
Weighted load-imbalance factor =   1.40%
---------------------------------------------------------------------------------------------------
```

Table format:
* `Time`: physical time
* `Step`: root-level updates
* `NPatch`: total number of patches
* `NPar`: total number of particles
* `Rank`: MPI rank
* `Level`: AMR level
* `X(Y)` in each row: `X` and `Y` are the estimated workload and load imbalance, respectively
* `Sum`: total workload on each level
* `Ave`: average workload on each level
* `Max`: maximum workload per rank on each level
* `Imb`: load imbalance on each level defined as `(Max-Ave)/Ave`
* `Weighted load-imbalance factor`: overall load imbalance, where the fact that different
levels may have different time-steps has been included in the weighting



<br>

## Links
* [[Simulation Logs]]