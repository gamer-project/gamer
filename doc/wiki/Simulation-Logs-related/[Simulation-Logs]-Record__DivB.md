This file records the numerical errors of the divergence-free constraint
on the magnetic field. The relative error is defined as
&Delta;&xi;<sub>l</sub>|&#8711;&#8901;<var>B</var>|/<|<var>B</var>|>,
where <var>B</var> is the face-centered magnetic field of a cell and
&Delta;&xi;<sub>l</sub> is the cell width along &xi; on level <var>l</var>.

Example:
``` markdown
# Tolerated Error: 9.9999997e-06

#-------------------------------------------------------------------------
#        Time         Step        AveError        MaxError     FailedCells
0.0000000e+00            0   0.0000000e+00   0.0000000e+00               0
1.7641872e-03            1   1.6332946e-08   5.8106949e-08               0
3.5260132e-03            2   2.1671077e-08   1.0531079e-07               0
5.2856565e-03            3   2.6458105e-08   1.1240161e-07               0
7.0427784e-03            4   3.0183041e-08   1.5688254e-07               0
8.7978840e-03            5   3.3473810e-08   1.6967800e-07               0
1.0550219e-02            6   3.6699682e-08   1.6826570e-07               0
1.2300350e-02            7   3.9571327e-08   2.0603242e-07               0
1.4047826e-02            8   4.2254163e-08   2.0667694e-07               0
1.5793136e-02            9   4.4946636e-08   2.2079354e-07               0
1.7535893e-02           10   4.7368575e-08   2.2033124e-07               0

```

Table format:
* `Time`: physical time
* `Step`: number of root-level updates
* `AveError`: average (L2-norm) error of all leaf cells
* `MaxError`: maximum error among all leaf cells
* `FailedCells`: number of cells exceeding the _Tolerated Error_
specified on the top of this file


> [!CAUTION]
> Errors should be on the order of machine precision but
may accumulate with time.

<br>

## Links
* [[Simulation Logs]]
