Parameters described on this page:
[OMP_NTHREAD](#OMP_NTHREAD), &nbsp;
[OPT__INIT_GRID_WITH_OMP](#OPT__INIT_GRID_WITH_OMP), &nbsp;
[LB_INPUT__WLI_MAX](#LB_INPUT__WLI_MAX), &nbsp;
[LB_INPUT__PAR_WEIGHT](#LB_INPUT__PAR_WEIGHT), &nbsp;
[OPT__RECORD_LOAD_BALANCE](#OPT__RECORD_LOAD_BALANCE), &nbsp;
[OPT__MINIMIZE_MPI_BARRIER](#OPT__MINIMIZE_MPI_BARRIER) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="OMP_NTHREAD"></a>
* #### `OMP_NTHREAD` &ensp; (&#8805;1; <1 &#8594; set to default) &ensp; [-1]
    * **Description:**
Number of OpenMP threads associated with each MPI process.
When enabling MPI, the default is set to the ratio between the total number of CPU cores
and the total number of MPI processes. When disabling MPI, the default
is set to the maximum number of threads available, which can be controlled
by the environment variable `OMP_NUM_THREADS`.
Having `OMP_NTHREAD=1` is equivalent to disabling the OpenMP parallelization.
See also [[Hybrid MPI/OpenMP/GPU | MPI-and-OpenMP#hybrid-mpiopenmpgpu]] and
[[MPI Binding and Thread Affinity | MPI-and-OpenMP#mpi-binding-and-thread-affinity]].
    * **Restriction:**
Only applicable when enabling the compilation option
[[--openmp | Installation:-Option-List#--openmp]].

<a name="OPT__INIT_GRID_WITH_OMP"></a>
* #### `OPT__INIT_GRID_WITH_OMP` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Whether or not to enable OpenMP when assigning the initial condition
of different grid patches. In can be enabled in most cases unless,
for example, the initial condition setup involves random numbers.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--openmp | Installation:-Option-List#--openmp]].

<a name="LB_INPUT__WLI_MAX"></a>
* #### `LB_INPUT__WLI_MAX` &ensp; (&#8805;0.0) &ensp; [0.1]
    * **Description:**
Weighted load imbalancing (WLI) threshold. Patches on all levels will
be redistributed among different MPI processes when the WLI factor is
estimated to be higher than a given threshold. See
[[Performance Optimizations: Load Balancing | Performance Optimizations:-Load-Balancing]]
for details.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--mpi | Installation:-Option-List#--mpi]].

<a name="LB_INPUT__PAR_WEIGHT"></a>
* #### `LB_INPUT__PAR_WEIGHT` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Load balancing weight of one particle over one cell. It is used
to improve load balancing for the simulations with particles. See
[[Performance Optimizations: Load Balancing | Performance Optimizations:-Load-Balancing]]
for details. The typical values are 1.0 ~ 2.0.
    * **Restriction:**
Only applicable when enabling the compilation options
[[--mpi | Installation:-Option-List#--mpi]] and
[[--particle | Installation:-Option-List#--particle]].

<a name="OPT__RECORD_LOAD_BALANCE"></a>
* #### `OPT__RECORD_LOAD_BALANCE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Record the load balancing information in the log file
[[Record__LoadBalance | Simulation-Logs:-Record__LoadBalance]].
    * **Restriction:**
Only applicable when enabling the compilation option
[[--mpi | Installation:-Option-List#--mpi]].

<a name="OPT__MINIMIZE_MPI_BARRIER"></a>
* #### `OPT__MINIMIZE_MPI_BARRIER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Minimize the MPI synchronization between grid and particle routines
to improve load balancing. It can improve the performance notably,
especially for simulations with particles.
    * **Restriction:**
For simulations with particles, one must enable the compilation option
[[--store_pot_ghost | Installation:-Option-List#--store_pot_ghost]] and
set [[PAR_IMPROVE_ACC | Runtime Parameters:-Particles#PAR_IMPROVE_ACC]]=1.
[[OPT__TIMING_BALANCE | Runtime Parameters:-Miscellaneous#OPT__TIMING_BALANCE]]
must be disabled. In addition, it is currently recommended to disable
[[AUTO_REDUCE_DT | Runtime Parameters:-Timestep#AUTO_REDUCE_DT]].


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of MPI and OpenMP | MPI-and-OpenMP]]
