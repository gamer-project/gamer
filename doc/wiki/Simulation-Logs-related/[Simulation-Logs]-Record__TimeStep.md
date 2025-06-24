This file records the timesteps set by various constraints.

Example:
```markdown
# Lv      Step   Counter        TimeOld        TimeNew          dTime      Hydro_CFL      Hydro_Acc      Data_Dump       End_Time        Par_Vel        Par_Acc      Sync_FaLv     Sync_SonLv      AutoRedDt
   0         0         0  0.0000000e+00  5.7810403e-02  5.7810403e-02  5.7810403e-02  3.8746542e-01  5.0000000e-01  1.0000000e+00            inf            inf  3.4028235e+38  0.0000000e+00  1.0000000e+00
   1         0         0  0.0000000e+00  2.7251148e-02  2.7251148e-02  2.8751362e-02  2.7327576e-01  5.0000000e-01  1.0000000e+00  2.7251148e-02  3.2869031e-01  5.7810403e-02  0.0000000e+00  1.0000000e+00
   2         0         0  0.0000000e+00  1.1388617e-02  1.1388617e-02  1.4373764e-02  1.8945998e-01  5.0000000e-01  1.0000000e+00  1.1388617e-02  1.3525940e-01  2.7251148e-02  3.4028235e+38  1.0000000e+00
   2         0         1  1.1388617e-02  2.2806077e-02  1.1417460e-02  1.4375145e-02  1.8895561e-01  4.8861138e-01  9.8861138e-01  1.1417460e-02  1.3519738e-01  1.5862531e-02  3.4028235e+38  1.0000000e+00
   2         0         2  2.2806077e-02  2.7251148e-02  4.4450713e-03  1.4371166e-02  1.8914053e-01  4.7719392e-01  9.7719392e-01  1.1404442e-02  1.3502470e-01  4.4450713e-03  3.4028235e+38  1.0000000e+00
   1         0         1  2.7251148e-02  5.4564332e-02  2.7313184e-02  2.8762037e-02  2.7449206e-01  4.7274885e-01  9.7274885e-01  2.7313184e-02  3.2819969e-01  3.0559255e-02  2.2808885e-02  1.0000000e+00
   2         0         3  2.7251148e-02  3.8643738e-02  1.1392589e-02  1.4369172e-02  1.8924589e-01  4.7274885e-01  9.7274885e-01  1.1392589e-02  1.3499083e-01  2.7313184e-02  3.4028235e+38  1.0000000e+00
   2         0         4  3.8643738e-02  5.0006536e-02  1.1362798e-02  1.4362061e-02  1.8951729e-01  4.6135626e-01  9.6135626e-01  1.1362798e-02  1.3511053e-01  1.5920594e-02  3.4028235e+38  1.0000000e+00
   2         0         5  5.0006536e-02  5.4564332e-02  4.5577961e-03  1.4361952e-02  1.8928616e-01  4.4999346e-01  9.4999346e-01  1.1333992e-02  1.3501261e-01  4.5577961e-03  3.4028235e+38  1.0000000e+00
   1         0         2  5.4564332e-02  5.7810403e-02  3.2460714e-03  2.8789273e-02  2.7419209e-01  4.4543567e-01  9.4543567e-01  2.7375017e-02  3.0900203e-01  3.2460714e-03  3.4028235e+38  1.0000000e+00
   2         0         6  5.4564332e-02  5.7810403e-02  3.2460714e-03  1.4362757e-02  1.8928155e-01  4.4543567e-01  9.4543567e-01  1.1322705e-02  1.3496141e-01  3.2460714e-03  3.4028235e+38  1.0000000e+00
```

Table format:
* `Lv`: AMR level
* `Step`: cumulative root-level updates
* `Counter`: cumulative updates on each level
* `TimeOld`, `TimeNew`: update this level from `TimeOld` to `TimeNew`=`TimeOld`+`dTime`
* `dTime`: adopted timestep (which is the minimum timestep of all constraints)
* `Hydro_CFL`: CFL condition in hydrodynamics (see
[[DT__FLUID | Runtime-Parameters:-Timestep#DT__FLUID]] and
[[DT__FLUID_INIT | Runtime-Parameters:-Timestep#DT__FLUID_INIT]])
* `Hydro_Acc`: gravitational acceleration in hydrodynamics (see
[[DT__GRAVITY | Runtime-Parameters:-Timestep#DT__GRAVITY]])
* `Data_Dump`: timestep to reach the next data dump time (see
[[OPT__OUTPUT_MODE | Runtime-Parameters:-Outputs#OPT__OUTPUT_MODE]] and
[[OUTPUT_DT | Runtime-Parameters:-Outputs#OUTPUT_DT]])
* `End_Time`: timestep to reach the simulation end time (see
[[END_T | Runtime-Parameters:-General#END_T]])
* `Par_Vel`: particle velocity (see
[[DT__PARVEL | Runtime-Parameters:-Timestep#DT__PARVEL]] and
[[DT__PARVEL_MAX | Runtime-Parameters:-Timestep#DT__PARVEL_MAX]])
* `Par_Acc`: particle acceleration (see
[[DT__PARACC | Runtime-Parameters:-Timestep#DT__PARACC]])
* `Sync_FaLv`: timestep to synchronize with the parent level (see also
[[DT__SYNC_PARENT_LV | Runtime-Parameters:-Timestep#DT__SYNC_PARENT_LV]])
* `Sync_SonLv`: timestep to help synchronize with the children level (see also
[[DT__SYNC_CHILDREN_LV | Runtime-Parameters:-Timestep#DT__SYNC_CHILDREN_LV]])
* `AutoRedDt`: timestep reducing factor used when the program fails (see
[[AUTO_REDUCE_DT | Runtime-Parameters:-Timestep#AUTO_REDUCE_DT]])

<br>

## Links
* [[Simulation Logs]]
