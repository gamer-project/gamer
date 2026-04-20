This file records the patch allocation and distribution.
The patch dimension is currently fixed to (8, 8, 8).

Example:
``` markdown
Time = 0.0000000e+00,  Step =       0,  NPatch =     258184

Rank        Level 0         Level 1         Level 2         Level 3         Level 4         Level 5
   0    928( 90.62%)   1088( 13.28%)   4376(  6.68%)  17360(  3.31%)  39160(  0.93%)      0(  0.00%)
   1   1120(109.38%)   1424( 17.38%)   5216(  7.96%)  18272(  3.49%)  40128(  0.96%)      0(  0.00%)
   2   1120(109.38%)   1432( 17.48%)   5224(  7.97%)  18272(  3.49%)  40136(  0.96%)      0(  0.00%)
   3    928( 90.62%)   1088( 13.28%)   4376(  6.68%)  17368(  3.31%)  39168(  0.93%)      0(  0.00%)
----------------------------------------------------------------------------------------------------
Sum:   4096(100.00%)   5032( 15.36%)  19192(  7.32%)  71272(  3.40%) 158592(  0.95%)      0(  0.00%)
Ave:   1024.00         1258.00         4798.00        17818.00        39648.00            0.00
Imb:   1120(  9.38%)   1432( 13.83%)   5224(  8.88%)  18272(  2.55%)  40136(  1.23%)      0(  0.00%)
Weighted load-imbalance factor =   1.70%
----------------------------------------------------------------------------------------------------

```

Table format:
* `Time`: physical time
* `Step`: number of root-level updates
* `NPatch`: total number of patches currently allocated on all levels in all MPI processes
* `Rank`: MPI rank
* `Level`: AMR level
* `X(Y)` in each row except for `Imb`: `X` is the number of patches and `Y` is the volume-filling fraction
* `Sum`: total number of patches on each level
* `Ave`: average number of patches per MPI process on each level
* `Imb`: each `X(Y)` shows the maximum number of patches in one MPI process (`X`)
and the corresponding load imbalance (`Y`), estimated by <var>(maximum-average)/average</var>
* `Weighted load-imbalance factor`: overall load imbalance, where different levels are weighted
by their number of updates (i.e., time-steps)

> [!CAUTION]
> Buffer patches (i.e., patches for filling the ghost zones) are not included here.

<br>

## Links
* [[Simulation Logs]]
