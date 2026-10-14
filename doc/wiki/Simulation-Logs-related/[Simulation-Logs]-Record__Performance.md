The file records the code performance.

Example:
```markdown
#         Time          Step               dt         NCell  NUpdate_Cell   ElapsedTime  Perf_Overall  Perf_PerRank     NParticle   NUpdate_Par  ParPerf_Overall  ParPerf_PerRank   NUpdate_Lv0   NUpdate_Lv1   NUpdate_Lv2   NUpdate_Lv3   NUpdate_Lv4   NUpdate_Lv5
 2.9144818e-03             1    2.9144818e-03      1.32e+08      2.29e+09      3.53e+01      6.48e+07      1.62e+07      6.68e+06      1.43e+08         4.04e+06         1.01e+06             1             2             4             8            24             0
 5.8224630e-03             2    2.9079812e-03      1.32e+08      2.30e+09      3.61e+01      6.36e+07      1.59e+07      6.68e+06      1.43e+08         3.95e+06         9.89e+05             1             2             4             8            24             0
 8.7230441e-03             3    2.9005811e-03      1.32e+08      2.29e+09      3.64e+01      6.28e+07      1.57e+07      6.68e+06      1.43e+08         3.92e+06         9.80e+05             1             2             4             8            24             0
 1.0000000e-02             4    1.2769559e-03      1.32e+08      1.06e+09      2.51e+01      4.24e+07      1.06e+07      6.68e+06      6.59e+07         2.62e+06         6.56e+05             1             1             2             4            11             0
```

Table format:
* `Time`: physical time
* `Step`: cumulative number of root-level updates
* `dt`: root-level time-step
* `NCell`: total number of cells
* `NUpdate_Cell`: total number of cell updates
* `ElapsedTime`: wall-clock time used by this step (in second)
* `Perf_Overall`: overall performance in cell updates per second
* `Perf_PerRank`: average performance per MPI process in cell updates per second
* `NParticle`: total number of particles
* `NUpdate_Par`: total number of particle updates
* `ParPerf_Overall`: overall performance in particle updates per second
* `ParPerf_PerRank`: average performance per MPI process in particle updates per second
* `NUpdate_Lv*`: number of time-steps on this level in this step

> [!NOTE]
> On Tesla K20 or K40 GPUs, the typical performance per GPU is 1e7 ~ 2e7
cell updates per second. Tesla P100 can be a factor of 2 faster.

<br>

## Links
* [[Simulation Logs]]
