Problem-specific runtime parameters. This file follows the same
syntax as [[Input__Parameter | Runtime-Parameters:-Input__Parameter#file-format]]. See
`example/test_problem/*/*/Input__TestProb` for some examples.
The following shows
`example/test_problem/Hydro/Riemann/Input__TestProb`:

```
Riemann_Prob     0     # target Riemann problem (0 : Sod's shock tube
                       #                         1 : strong shock
                       #                         2 : two shocks
                       #                         3 : Einfeldt's 1-2-0-3
                       #                         4 : Einfeldt's 1-1-2-5
                       #                         5 : sonic rarefaction wave
Riemann_LR      +1     # wave propagation direction (>0/<0 --> positive/negative direction) [1]
Riemann_XYZ      0     # wave propagation axis (0/1/2 --> x/y/z) [0]
```
