## Running in Parallel

### MPI-only
To enable MPI support, follow the steps below:
1. Edit the `Makefile` and recompile the code (see [[Installation]] for details)
    1. Set `CXX` to an MPI compiler (e.g., `mpicxx`)
    2. Set `MPI_PATH` to you MPI installation path
    3. Turn off
[[SERIAL | Installation: Simulation-Options#SERIAL]]
and turn on
[[LOAD_BALANCE=HILBERT | Installation: Simulation-Options#LOAD_BALANCE]]
    4. Recompile the code by `make clean; make`

2. Launch the code with MPI (consult your system documentation),
for instance,

    ```bash
    mpirun -np 10 ./gamer
    ```

### Hybrid MPI/OpenMP
To enable hybrid MPI/OpenMP, follow the MPI-only prescriptions
given above with the following additional steps:
1. Also turn on the compilation option
[[OPENMP | Installation: Simulation-Options#OPENMP]]
and set the OpenMP flag `OPENMPFLAG` properly in the `Makefile`

2. Set the number of threads through the runtime parameter
[[OMP_NTHREAD | MPI-and-OpenMP#OMP_NTHREAD]]

3. The recommended way for launching hybrid MPI/OpenMP jobs can vary from
system to system. So please consult your system documentation. Also, check
out the [Remarks](#remarks) below.


## Compilation Options

Related options:
[[SERIAL | Installation: Simulation-Options#SERIAL]], &nbsp;
[[LOAD_BALANCE | Installation: Simulation-Options#LOAD_BALANCE]], &nbsp;
[[OPENMP | Installation: Simulation-Options#OPENMP]] &nbsp;


## Runtime Parameters

Parameters described on this page:
[OMP_NTHREAD](#OMP_NTHREAD), &nbsp;
[OPT__INIT_GRID_WITH_OMP](#OPT__INIT_GRID_WITH_OMP), &nbsp;
[LB_INPUT__WLI_MAX](#LB_INPUT__WLI_MAX), &nbsp;
[LB_INPUT__PAR_WEIGHT](#LB_INPUT__PAR_WEIGHT), &nbsp;
[OPT__RECORD_LOAD_BALANCE](#OPT__RECORD_LOAD_BALANCE), &nbsp;
[OPT__MINIMIZE_MPI_BARRIER](#OPT__MINIMIZE_MPI_BARRIER) &nbsp;

Other related parameters: none

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
See also [Hybrid MPI/OpenMP/GPU](#hybrid-mpiopenmpgpu) and
[MPI Binding and Thread Affinity](#mpi-binding-and-thread-affinity).
    * **Restriction:**
Only applicable when enabling the compilation option
[[OPENMP | Installation: Simulation-Options#OPENMP]].

<a name="OPT__INIT_GRID_WITH_OMP"></a>
* #### `OPT__INIT_GRID_WITH_OMP` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Whether or not to enable OpenMP when assigning the initial condition
of different grid patches. In can be enabled in most cases unless,
for example, the initial condition setup involves random numbers.
    * **Restriction:**
Only applicable when enabling the compilation option
[[OPENMP | Installation: Simulation-Options#OPENMP]].

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
[[LOAD_BALANCE | Installation: Simulation-Options#LOAD_BALANCE]].

<a name="LB_INPUT__PAR_WEIGHT"></a>
* #### `LB_INPUT__PAR_WEIGHT` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Load balancing weight of one particle over one cell. It is used
to improve load balancing for the simulations with particles. See
[[Performance Optimizations: Load Balancing | Performance Optimizations:-Load-Balancing]]
for details. The typical values are 1.0 ~ 2.0.
    * **Restriction:**
Only applicable when enabling the compilation options
[[LOAD_BALANCE | Installation: Simulation-Options#LOAD_BALANCE]] and
[[PARTICLE | Installation: Simulation-Options#PARTICLE]].

<a name="OPT__RECORD_LOAD_BALANCE"></a>
* #### `OPT__RECORD_LOAD_BALANCE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Record the load balancing information in the log file
[[Record__LoadBalance | Simulation-Logs:-Record__LoadBalance]].
    * **Restriction:**
Only applicable when enabling the compilation option
[[LOAD_BALANCE | Installation: Simulation-Options#LOAD_BALANCE]].

<a name="OPT__MINIMIZE_MPI_BARRIER"></a>
* #### `OPT__MINIMIZE_MPI_BARRIER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Minimize the MPI synchronization between grid and particle routines
to improve load balancing. It can improve the performance notably,
especially for simulations with particles.
    * **Restriction:**
For simulations with particles, one must enable the compilation option
[[STORE_POT_GHOST | Installation: Simulation-Options#STORE_POT_GHOST]] and
set [[PAR_IMPROVE_ACC | Particles#PAR_IMPROVE_ACC]]=1.
[[OPT__TIMING_BALANCE | Runtime Parameters:-Miscellaneous#OPT__TIMING_BALANCE]]
must be disabled. In addition, it is currently recommended to disable
[[AUTO_REDUCE_DT | Runtime Parameters:-Timestep#AUTO_REDUCE_DT]].


## Remarks

### Hybrid MPI/OpenMP/GPU

The code does not pose any restrictions on the ratio between the number
of MPI processes per node and the number of OpenMP threads associated
with each MPI process. When using GPUs, the most straightforward way is
to set the total number of MPI processes equal
to the total number of GPUs you want to use, and then set the number of
OpenMP threads per process ([OMP_NTHREAD](#OMP_NTHREAD)) equal to the
ratio between the number of CPU cores per node and the number of GPUs per
node.

For example, to run a job using 8 nodes, each of which is composed of 2
ten-core CPUs and 2 GPUs, you could consider adopting 16 MPI processes and
10 OpenMP threads. This will allow different MPI processes to use different
ten-core CPUs and GPUs, and different OpenMP threads of the same MPI process
to use different CPU cores of the same CPU.

Depending on the system specifications, it might be beneficial to allow
multiple MPI processes to access the same GPU. Using CUDA
[Multi-Process Service (MPS)](https://docs.nvidia.com/deploy/pdf/CUDA_Multi_Process_Service_Overview.pdf)
is recommended for this case and may boost the performance further.
For instance, for the example above, you could adopt 32 MPI processes and
5 OpenMP threads. This will allow two MPI processes to share the same GPUs
while still allow different OpenMP threads of the same MPI process to use
different CPU cores of the same CPU.
See also [[Set and Validate GPU IDs | GPU#set-and-validate-gpu-ids]].

Please also check
[MPI Binding and Thread Affinity](#mpi-binding-and-thread-affinity)
carefully.


### MPI Binding and Thread Affinity

It is important to ensure that different MPI processes and OpenMP
threads access different CPU cores. The recommended setup can vary
from system to system. So please check your system documentation.
One straightforward way to validate it is to use the Linux command
`top` and type `1` to check whether the CPU usage of ALL cores is
close to 100%.

One can also validate the MPI and OpenMP binding by searching for the
keyword "OpenMP" in the log file `Record__Note`. The following example
adopts 8 MPI processes and `OMP_NTHREAD=10` to run a job on 4 nodes
named golub121-124, each of which is composed of 2 ten-core CPUs and 2 GPUs:
<pre>
OpenMP Diagnosis
***********************************************************************************
OMP__SCHEDULE                   DYNAMIC
OMP__SCHEDULE_CHUNK_SIZE        1
OMP__NESTED                     OFF

CPU core IDs of all OpenMP threads (tid == thread ID):
------------------------------------------------------------------------
 Rank        Host  NThread  tid-00  tid-01  tid-02  tid-03  tid-04  tid-05  tid-06  tid-07  tid-08  tid-09
    0    golub121       10       0       2       4       6       8      10      12      14      16      18
    1    golub121       10       1       3       5       7       9      11      13      15      17      19
    2    golub122       10       0       2       4       6       8      10      12      14      16      18
    3    golub122       10       1       3       5       7       9      11      13      15      17      19
    4    golub123       10       0       2       4       6       8      10      12      14      16      18
    5    golub123       10       1       3       5       7       9      11      13      15      17      19
    6    golub124       10       0       2       4       6       8      10      12      14      16      18
    7    golub124       10       1       3       5       7       9      11      13      15      17      19
***********************************************************************************
</pre>
Check the following things:
* **The number under `NThread` is the same as the runtime parameter
[OMP_NTHREAD](#OMP_NTHREAD)**
* **Different OpenMP threads use different CPU cores**

To achieve an optimal performance, it is also important to take into
account thread affinity and non-uniform memory access (NUMA).
Generally, it is recommended to have OpenMP threads running in the
same NUMA domain to improve memory affinity. But one still needs
to experiment with different configurations to fine-tune the performance.
The Linux command `lscpu` can be used to display information about
your CPU architecture. For example, on a node with 2 ten-core CPUs
(as the example given above), it shows
<pre>
...
CPU(s):                20
On-line CPU(s) list:   0-19
Thread(s) per core:    1
Core(s) per socket:    10
Socket(s):             2
NUMA node(s):          2
...
NUMA node0 CPU(s):     0,2,4,6,8,10,12,14,16,18
NUMA node1 CPU(s):     1,3,5,7,9,11,13,15,17,19
</pre>
By comparing it with the thread binding information recorded in
the log file `Record__Note` (as described above), it confirms that
in this example different threads in the same MPI process do run
in the same NUMA domain.


### OpenMP Support in GRACKLE
See [[Library Configurations -- GRACKLE | Installation:-External-Libraries#grackle]]
for how to enable OpenMP in GRACKLE.


<br>

## Links
* [[How to run the code | Running the code]]
* [[How to install GRACKLE | Installation:-External-Libraries#grackle]]
* [[Main page of Runtime Parameters | Runtime Parameters]]
