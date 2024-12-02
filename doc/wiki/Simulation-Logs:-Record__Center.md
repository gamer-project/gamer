This file records the position of maximum density, minimum potential, and center of mass at each step.

Example:
``` markdown
#               Time        Step         MaxDens       MaxDens_x       MaxDens_y       MaxDens_z      MaxParDens    MaxParDens_x    MaxParDens_y    MaxParDens_z    MaxTotalDens  MaxTotalDens_x  MaxTotalDens_y  MaxTotalDens_z         MinPote       MinPote_x       MinPote_y       MinPote_z     Final_NIter        Final_dR           CoM_x           CoM_y           CoM_z
0.00000000000000e+00           0   2.4367696e-01   1.4941406e+00   1.4941406e+00   1.4941406e+00   9.1377342e-01   1.5058594e+00   1.5058594e+00   1.5058594e+00   1.1574503e+00   1.5058594e+00   1.5058594e+00   1.5058594e+00  -4.1190378e-02   1.5058594e+00   1.4941406e+00   1.5058594e+00               1   9.6799215e-03   1.4995081e+00   1.5005405e+00   1.5008522e+00
9.00683318704315e-02           1   2.4402954e-01   1.5058594e+00   1.5058594e+00   1.5058594e+00   8.8797969e-01   1.5175781e+00   1.5058594e+00   1.5058594e+00   1.1314732e+00   1.5058594e+00   1.5058594e+00   1.5058594e+00  -4.1192625e-02   1.5058594e+00   1.5058594e+00   1.5058594e+00               1   7.5696964e-03   1.5013204e+00   1.5023092e+00   1.5009508e+00
1.81193407302794e-01           2   2.4525556e-01   1.5058594e+00   1.5058594e+00   1.5058594e+00   8.5094547e-01   1.5058594e+00   1.5058594e+00   1.4824219e+00   1.0803783e+00   1.5058594e+00   1.5058594e+00   1.4824219e+00  -4.1127916e-02   1.5058594e+00   1.5058594e+00   1.5058594e+00               1   1.8903390e-02   1.5031526e+00   1.5041017e+00   1.5010477e+00
2.73344760914069e-01           3   2.4442711e-01   1.5058594e+00   1.5058594e+00   1.5058594e+00   8.5942018e-01   1.5058594e+00   1.4941406e+00   1.4941406e+00   1.0986829e+00   1.5058594e+00   1.4941406e+00   1.4941406e+00  -4.1115854e-02   1.5058594e+00   1.5058594e+00   1.5058594e+00               1   1.3729049e-02   1.5050049e+00   1.5059173e+00   1.5011455e+00
3.65979077563538e-01           4   2.4281065e-01   1.5058594e+00   1.5058594e+00   1.4941406e+00   8.4279549e-01   1.4941406e+00   1.4941406e+00   1.5058594e+00   1.0680594e+00   1.4941406e+00   1.4941406e+00   1.5058594e+00  -4.1099463e-02   1.5058594e+00   1.5058594e+00   1.4941406e+00               1   1.9192376e-02   1.5068745e+00   1.5077393e+00   1.5012472e+00
4.56276972711556e-01           5   2.3963770e-01   1.5058594e+00   1.5058594e+00   1.4941406e+00   8.4634453e-01   1.5175781e+00   1.4941406e+00   1.4824219e+00   1.0625260e+00   1.5175781e+00   1.4941406e+00   1.4824219e+00  -4.1039847e-02   1.5058594e+00   1.5058594e+00   1.4941406e+00               1   2.5946765e-02   1.5086964e+00   1.5095107e+00   1.5013457e+00
```

Table format:
* `Time`: physical time
* `Step`: cumulative root-level updates
* `MaxDens`: value of the maximum fluid density
* `MaxDens_x/y/z`: x/y/z coordinate of the maximum fluid density position
* `MaxParDens`: value of the maximum particle density
* `MaxParDens_x/y/z`: x/y/z coordinate of the maximum particle density position
* `MaxTotalDens`: value of the maximum fluid+particle density
* `MaxTotalDens_x/y/z`: x/y/z coordinate of the maximum fluid+particle density position
* `MinPote`: value of the minimum gravitational potential
* `MinPote_x/y/z`: x/y/z coordinate of the minimum gravitational potential position
* `Final_NIter`: total number of iterations for finding the center of mass (see also
[[COM_MAX_ITER | Runtime-Parameters:-Miscellaneous#COM_MAX_ITER]])
* `Final_dR`: distance of the center coordinates update in the last iteration for finding the center of mass (see also
[[COM_TOLERR_R | Runtime-Parameters:-Miscellaneous#COM_TOLERR_R]])
* `CoM_x/y/z`: x/y/z coordinate of the center of mass (see also
[[COM_CEN_X | Runtime-Parameters:-Miscellaneous#COM_CEN_X]],
[[COM_CEN_Y | Runtime-Parameters:-Miscellaneous#COM_CEN_Y]],
[[COM_CEN_Z | Runtime-Parameters:-Miscellaneous#COM_CEN_Z]],
[[COM_MAX_R | Runtime-Parameters:-Miscellaneous#COM_MAX_R]],
[[COM_MIN_RHO | Runtime-Parameters:-Miscellaneous#COM_MIN_RHO]])


<br>

## Links
* [[Simulation Logs]]