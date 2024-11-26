This file records the CPU memory consumption.

Example:
``` markdown
#         Time          Step     Virtual_Max (MB)    Virtual_Sum (MB)   Resident_Max (MB)   Resident_Sum (MB)
 0.0000000e+00             0           102868.36           205735.19             1741.32             3481.74
 5.7810403e-02             1           102868.57           205735.47             1749.68             3498.42
 1.1256044e-01             2           102868.61           205735.53             1749.91             3498.71
 1.7039426e-01             3           102873.61           205742.53             1754.65             3504.28
 2.2822480e-01             4           102873.61           205742.53             1754.67             3504.43
 2.7323392e-01             5           102876.99           205746.29             1757.69             3507.81
```

Table format:
* `Time`: physical time
* `Step`: number of root-level updates
* `Virtual_Max`: maximum virtual memory consumption in one MPI process
* `Virtual_Sum`: total virtual memory consumption in all MPI processes
* `Resident_Max`: maximum resident memory consumption in one MPI process
* `Resident_Sum`: total resident memory consumption in all MPI processes


> [!CAUTION]
> CUDA may consume very large virtual memory. However, it doesn't seem to be
a serious problem in most cases.

<br>

## Links
* [[Simulation Logs]]
