This page describes various data formats and options for storing simulation snapshots.
For simulation logs recording, for example, performance, load balancing, and memory
consumption, see [[Simulation Logs | Simulation-Logs]].


## Compilation Options

Related options:
[[--hdf5 | [Installation]-Option-List#--hdf5]] &nbsp;


## Runtime Parameters
[[ [Runtime parameters] Outputs | [Runtime-Parameters]-Outputs ]]

Other related parameters: none


## Remarks

### Enforcing Outputs

Typically, the simulation output interval is controlled by
[[OUTPUT_STEP | [Runtime-Parameters]-Outputs#OUTPUT_STEP]],
[[OUTPUT_DT | [Runtime-Parameters]-Outputs#OUTPUT_DT]], or the timetable
[[Input__DumpTable | [Runtime-Parameters]-Input__DumpTable]]
depending on
[[OPT__OUTPUT_MODE | [Runtime-Parameters]-Outputs#OPT__OUTPUT_MODE]].
One can also enforce output by creating a file named
`DUMP_GAMER_DUMP` (e.g., using `touch DUMP_GAMER_DUMP`) in the same
directory as the executable `gamer`. The program will check the existence
of this file after every root-level update. If it is detected,
the program will output data (assuming at least one of the output options described
on this page is enabled) and delete the file `DUMP_GAMER_DUMP`.
Note that this functionality is supported only when the runtime parameter
[[OPT__MANUAL_CONTROL | [Runtime-Parameters]-Miscellaneous#OPT__MANUAL_CONTROL]]
is on.


### Sub-cadence Outputs

The runtime parameter
[[OPT__OUTPUT_SUBDIV | [Runtime-Parameters]-Outputs#OPT__OUTPUT_SUBDIV]]
fires additional outputs at a finer cadence than (and between) the main data dumps
without writing extra `Data_XXXXXX` files. Sub-cadence data are written to a single
`SubData_*` file in the full-snapshot HDF5 format, whose `Tree`, `GridData` (field list
set by `Input__Sub_Grid`), and `Particle` (massive and/or tracer particles, with optional
float32 output) groups are controlled independently, from fully `yt`-loadable snapshots
down to minimal particles-only files for direct `h5py` access; sub-cadence
`Output_User_Ptr()` calls are also supported. All sub-cadence outputs are indexed by the
global counter `SubDumpID` and recorded in
[[Record__TimeSubDump | [Simulation-Logs]-Record__TimeSubDump]].


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Simulation Logs | Simulation-Logs]]
