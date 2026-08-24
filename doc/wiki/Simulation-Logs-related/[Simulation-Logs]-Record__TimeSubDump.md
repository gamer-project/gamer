This file records the physical time and root-level steps associated with
each sub-cadence output (see
[[OPT__OUTPUT_SUBDIV | [Runtime-Parameters]-Outputs#OPT__OUTPUT_SUBDIV]]),
including the sub-dump events coinciding with the main data dumps.

Example:
``` markdown
SubDumpID                       Time                 Step
        0       0.00000000000000e+00                    0
        1       2.50000000000000e-01                    5
        2       5.00000000000000e-01                   10
```

Table format:
* `SubDumpID`: sub-dump ID (distinct from `DumpID`; advances at every sub-dump event)
* `Time`: physical time
* `Step`: number of root-level updates


<br>

## Links
* [[Simulation Logs]]
