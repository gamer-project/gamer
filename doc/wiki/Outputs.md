This page describes various data formats and options for storing simulation snapshots.
For simulation logs recording, for example, performance, load balancing, and memory
consumption, see [[Simulation Logs | Simulation-Logs]].


## Compilation Options

Related options:
[[--hdf5 | Installation:-Option-List#--hdf5]] &nbsp;


## Runtime Parameters
[[Runtime parameters: Outputs | Runtime-Parameters:-Outputs]]

Other related parameters: none


## Remarks

### Enforcing Outputs

Typically, the simulation output interval is controlled by
[[OUTPUT_STEP | Runtime-Parameters:-Outputs#OUTPUT_STEP]],
[[OUTPUT_DT | Runtime-Parameters:-Outputs#OUTPUT_DT]], or the timetable
[[Input__DumpTable | Runtime-Parameters:-Input__DumpTable]]
depending on
[[OPT__OUTPUT_MODE | Runtime-Parameters:-Outputs#OPT__OUTPUT_MODE]].
One can also enforce output by creating a file named
`DUMP_GAMER_DUMP` (e.g., using `touch DUMP_GAMER_DUMP`) in the same
directory as the executable `gamer`. The program will check the existence
of this file after every root-level update. If it is detected,
the program will output data (assuming at least one of the output options described
on this page is enabled) and delete the file `DUMP_GAMER_DUMP`.
Note that this functionality is supported only when the runtime parameter
[[OPT__MANUAL_CONTROL | Runtime-Parameters:-Miscellaneous#OPT__MANUAL_CONTROL]]
is on.


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Simulation Logs | Simulation-Logs]]
