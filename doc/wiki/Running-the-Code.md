This page covers the following topics:
* [Running Simulations](running-simulations)
* [Restarting from a Snapshot](#restarting-from-a-snapshot)
* [Terminating or Pausing Simulations](#terminating-or-pausing-simulations)
* [Taking Notes](#taking-notes)

It is highly recommended that you first follow the demos given in
[[Quick Start | Quick-Start]].


## Running Simulations

To run the code in a serial mode (i.e., no MPI and OpenMP),
follow the steps below:

1. Compile the code (See [[Installation]] for details)
    1. Set `CXX` in the [[configuration file | Installation:-Machine-Configuration-File]] to a serial compiler (e.g., `g++` or `icpc`)
    2. Set the following arguments when generating `Makefile`:
       * [[--mpi | Installation:-Option-List#--mpi]]=`false`
       * [[--openmp | Installation:-Option-List#--openmp]]=`false`
    3. Recompile the code by `make clean; make`

2. Set the runtime parameters (see [[Runtime Parameters]] for details)
    1. Put the required input files (e.g., `Input__Parameter`) in the same
directory as the executable `gamer`
    2. Set parameters properly

3. Launch the code

    ```bash
    ./gamer
    ```

4. Check the outputs
    * Simulation progress &#8594; stdout
    * Simulation [[log files | Outputs#log-files]] &#8594; `Record__*`
    * Simulation [[snapshots | Outputs#snapshots]] &#8594; `Data_??????`

Further readings:
* [[Enable MPI and OpenMP parallelization | MPI and OpenMP]]
* [[Enable GPU acceleration | GPU]]


## Restarting from a Snapshot

To restart a simulation from a snapshot (e.g., `Data_000123`),
do the following steps:
1. Set [[OPT__INIT | Runtime-Parameters:-Initial-Conditions#OPT__INIT]]=2
2. Create a soft link named `RESTART` to a targeted snapshot:
`ln -s Data_000123 RESTART`
3. Run the code as usual

> [!NOTE]
> Many runtime parameters can be changed during restart
(e.g., number of MPI processes and OpenMP threads, maximum refinement level,
refinement criteria, time-step criteria). The next data dump ID will be
automatically set to the ID of the restart file plus one (unless specified by
[[INIT_DUMPID | Runtime-Parameters:-Outputs#INIT_DUMPID]]).


## Terminating or Pausing Simulations

Typically, a simulation ends when reaching either
[[END_T | Runtime-Parameters:-General#END_T]]
or
[[END_STEP | Runtime-Parameters:-General#END_STEP]].
One can also forcibly terminate a simulation by creating a file named
`STOP_GAMER_STOP` (e.g., using `touch STOP_GAMER_STOP`) in the same
directory as the executable `gamer`. The program will check the existence
of this file after every root-level update. If it is detected, the
program will output the final data (assuming at least one of the output options
described in [[Outputs]] is enabled), delete the file `STOP_GAMER_STOP`,
and then be terminated. Note that this functionality is supported only
when the runtime parameter
[[OPT__MANUAL_CONTROL | Runtime-Parameters:-Miscellaneous#OPT__MANUAL_CONTROL]]
is on.

Similarly, one can create a file named `PAUSE_GAMER_PAUSE` to pause a simulation.
To resume it, simply delete this file.


## Taking Notes

One can put simulation notes in a text file named `Input__Note`. The content
of this file will be automatically copied to the top of the log file
[[Record__Note | Simulation Logs:-Record__Note]] during code initialization. For example,
```bash
cat Input__Note
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
~ PUT YOUR SIMULATION NOTES HERE ~
</pre>
</details>

```bash
./gamer
head Record__Note
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
Simulation Note
***********************************************************************************

~ PUT YOUR SIMULATION NOTES HERE ~

***********************************************************************************
</pre>
</details>

<br>

## Links
* [[Installation]]
* [[Runtime Parameters]]
