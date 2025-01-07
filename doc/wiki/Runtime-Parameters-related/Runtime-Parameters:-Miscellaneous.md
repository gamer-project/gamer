Parameters described on this page:
[OPT__VERBOSE](#OPT__VERBOSE), &nbsp;
[OPT__TIMING_BARRIER](#OPT__TIMING_BARRIER), &nbsp;
[OPT__TIMING_BALANCE](#OPT__TIMING_BALANCE), &nbsp;
[OPT__TIMING_MPI](#OPT__TIMING_MPI), &nbsp;
[OPT__RECORD_NOTE](#OPT__RECORD_NOTE), &nbsp;
[OPT__RECORD_UNPHY](#OPT__RECORD_UNPHY), &nbsp;
[OPT__RECORD_MEMORY](#OPT__RECORD_MEMORY), &nbsp;
[OPT__RECORD_PERFORMANCE](#OPT__RECORD_PERFORMANCE), &nbsp;
[OPT__MANUAL_CONTROL](#OPT__MANUAL_CONTROL), &nbsp;
[OPT__RECORD_CENTER](#OPT__RECORD_CENTER), &nbsp;
[COM_CEN_X](#COM_CEN_X), &nbsp;
[COM_CEN_Y](#COM_CEN_Y), &nbsp;
[COM_CEN_Z](#COM_CEN_Z), &nbsp;
[COM_MAX_R](#COM_MAX_R), &nbsp;
[COM_MIN_RHO](#COM_MIN_RHO), &nbsp;
[COM_TOLERR_R](#COM_TOLERR_R), &nbsp;
[COM_MAX_ITER](#COM_MAX_ITER), &nbsp;
[OPT__RECORD_USER](#OPT__RECORD_USER), &nbsp;
[OPT__OPTIMIZE_AGGRESSIVE](#OPT__OPTIMIZE_AGGRESSIVE), &nbsp;
[OPT__SORT_PATCH_BY_LBIDX](#OPT__SORT_PATCH_BY_LBIDX), &nbsp;
[OPT__CK_REFINE](#OPT__CK_REFINE), &nbsp;
[OPT__CK_PROPER_NESTING](#OPT__CK_PROPER_NESTING), &nbsp;
[OPT__CK_CONSERVATION](#OPT__CK_CONSERVATION), &nbsp;
[OPT__CK_NORMALIZE_PASSIVE](#OPT__CK_NORMALIZE_PASSIVE), &nbsp;
[OPT__CK_RESTRICT](#OPT__CK_RESTRICT), &nbsp;
[OPT__CK_FINITE](#OPT__CK_FINITE), &nbsp;
[OPT__CK_PATCH_ALLOCATE](#OPT__CK_PATCH_ALLOCATE), &nbsp;
[OPT__CK_FLUX_ALLOCATE](#OPT__CK_FLUX_ALLOCATE), &nbsp;
[OPT__CK_NEGATIVE](#OPT__CK_NEGATIVE), &nbsp;
[OPT__CK_MEMFREE](#OPT__CK_MEMFREE), &nbsp;
[OPT__CK_PARTICLE](#OPT__CK_PARTICLE), &nbsp;
[OPT__CK_INTERFACE_B](#OPT__CK_INTERFACE_B), &nbsp;
[OPT__CK_DIVERGENCE_B](#OPT__CK_DIVERGENCE_B), &nbsp;
[OPT__CK_INPUT_FLUID](#OPT__CK_INPUT_FLUID) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="OPT__VERBOSE"></a>
* #### `OPT__VERBOSE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Show details of the simulation progress.
    * **Restriction:**

<a name="OPT__TIMING_BARRIER"></a>
* #### `OPT__TIMING_BARRIER` &ensp; (0=off, 1=on; <0 &#8594; set to default) &ensp; [-1]
    * **Description:**
Synchronize all MPI processes (by invoking `MPI_Barrier()`) before timing. It will lead to
more accurate timing results but may also deteriorate the performance.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--timing | Installation:-Option-List#--timing]].

<a name="OPT__TIMING_BALANCE"></a>
* #### `OPT__TIMING_BALANCE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Record the maximum and minimum times consumed by various major routines among
all MPI processes to check load balancing. The results will be recorded in the file
[[Record__Timing | Simulation-Logs:-Record__Timing]].
    * **Restriction:**
Only applicable when enabling the compilation option
[[--timing | Installation:-Option-List#--timing]].

<a name="OPT__TIMING_MPI"></a>
* #### `OPT__TIMING_MPI` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Record the MPI bandwidth achieved by various MPI calls in the file
[[Record__TimingMPI_Rank* | Simulation-Logs:-Record__TimingMPI_Rank*]].
    * **Restriction:**
Only applicable when enabling both the compilation options
[[--timing | Installation:-Option-List#--timing]]
and
[[--mpi | Installation:-Option-List#--mpi]].

<a name="OPT__RECORD_NOTE"></a>
* #### `OPT__RECORD_NOTE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Record the general simulation information, such as compilation options and
runtime parameters, in the file
[[Record__Note | Simulation-Logs:-Record__Note]].
    * **Restriction:**

<a name="OPT__RECORD_UNPHY"></a>
* #### `OPT__RECORD_UNPHY` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Record the number of cells with unphysical results being corrected in the file
[[Record__NCorrUnphy | Simulation-Logs:-Record__NCorrUnphy]].
    * **Restriction:**

<a name="OPT__RECORD_MEMORY"></a>
* #### `OPT__RECORD_MEMORY` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Record the total memory consumption of each MPI process in the file
[[Record__MemInfo | Simulation-Logs:-Record__MemInfo]].
    * **Restriction:**

<a name="OPT__RECORD_PERFORMANCE"></a>
* #### `OPT__RECORD_PERFORMANCE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Record the code performance in the file
[[Record__Performance | Simulation-Logs:-Record__Performance]].
    * **Restriction:**
Only applicable when enabling the compilation option
[[--timing | Installation:-Option-List#--timing]].

<a name="OPT__MANUAL_CONTROL"></a>
* #### `OPT__MANUAL_CONTROL` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Support enforcing data dump and/or terminating simulations during the runtime. See
[[ Enforcing Outputs | Outputs#enforcing-outputs]] and
[[Terminating Simulations | Running-the-Code#terminating-simulations]]
for details.
    * **Restriction:**

<a name="OPT__RECORD_CENTER"></a>
* #### `OPT__RECORD_CENTER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Record the position of maximum density, minimum potential, and center of mass in the file
[[Record__Center | Simulation-Logs:-Record__Center]].
    * **Restriction:**

<a name="COM_CEN_X"></a>
* #### `COM_CEN_X` &ensp; (within the simulation domain; if one of COM_CEN_X/Y/Z < 0.0 -> peak density position x) &ensp; [-1.0]
    * **Description:**
x coordinate as an initial guess for determining the center of mass for [OPT__RECORD_CENTER](#OPT__RECORD_CENTER).
    * **Restriction:**

<a name="COM_CEN_Y"></a>
* #### `COM_CEN_Y` &ensp; (within the simulation domain; if one of COM_CEN_X/Y/Z < 0.0 -> peak density position y) &ensp; [-1.0]
    * **Description:**
y coordinate as an initial guess for determining the center of mass for [OPT__RECORD_CENTER](#OPT__RECORD_CENTER).
    * **Restriction:**

<a name="COM_CEN_Z"></a>
* #### `COM_CEN_Z` &ensp; (within the simulation domain; if one of COM_CEN_X/Y/Z < 0.0 -> peak density position z) &ensp; [-1.0]
    * **Description:**
z coordinate as an initial guess for determining the center of mass for [OPT__RECORD_CENTER](#OPT__RECORD_CENTER).
    * **Restriction:**

<a name="COM_MAX_R"></a>
* #### `COM_MAX_R` &ensp; (>0.0; <0.0 &#8594; `__FLT_MAX__`) &ensp; [-1.0]
    * **Description:**
Maximum radius for determining the center of mass for [OPT__RECORD_CENTER](#OPT__RECORD_CENTER).
    * **Restriction:**

<a name="COM_MIN_RHO"></a>
* #### `COM_MIN_RHO` &ensp; (&#8805;0.0) &ensp; [-1.0]
    * **Description:**
Minimum density for determining the center of mass for [OPT__RECORD_CENTER](#OPT__RECORD_CENTER).
    * **Restriction:**

<a name="COM_TOLERR_R"></a>
* #### `COM_TOLERR_R` &ensp; (>0.0; <0.0 &#8594; maximum spatial resolution) &ensp; [-1.0]
    * **Description:**
Maximum tolerated error of deviation in distance during the iterations of determining the center of mass
for [OPT__RECORD_CENTER](#OPT__RECORD_CENTER).
    * **Restriction:**

<a name="COM_MAX_ITER"></a>
* #### `COM_MAX_ITER` &ensp; (&#8805;1) &ensp; [10]
    * **Description:**
Maximum number of iterations for determining the center of mass for [OPT__RECORD_CENTER](#OPT__RECORD_CENTER).
    * **Restriction:**

<a name="OPT__RECORD_USER"></a>
* #### `OPT__RECORD_USER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Perform user-specified simulation diagnostics. See
[[Problem Specific Diagnostics | Adding-New-Simulations#problem-specific-diagnostics]]
for details.
    * **Restriction:**

<a name="OPT__OPTIMIZE_AGGRESSIVE"></a>
* #### `OPT__OPTIMIZE_AGGRESSIVE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Apply aggressive performance optimizations. This option is experimental.
    * **Restriction:**

<a name="OPT__SORT_PATCH_BY_LBIDX"></a>
* #### `OPT__SORT_PATCH_BY_LBIDX` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Sort patches by their load-balance indices before dumping data and during restart
to improve bitwise reproducibility, especially when restarting with the same number
of MPI processes.
    * **Restriction:**
Not supported by [[--mpi | Installation:-Option-List#--mpi]]=false.

<a name="OPT__CK_REFINE"></a>
* #### `OPT__CK_REFINE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check the grid refinement (deprecated).
    * **Restriction:**
Only work with the refinement option
[[OPT__FLAG_RHO| Runtime-Parameters:-Refinement#OPT__FLAG_RHO]].

<a name="OPT__CK_PROPER_NESTING"></a>
* #### `OPT__CK_PROPER_NESTING` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check the proper-nesting condition.
    * **Restriction:**

<a name="OPT__CK_CONSERVATION"></a>
* #### `OPT__CK_CONSERVATION` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check the conservation laws. The results will be recorded in the file
[[Record__Conservation | Simulation-Logs:-Record__Conservation]].
    * **Restriction:**

<a name="OPT__CK_NORMALIZE_PASSIVE"></a>
* #### `OPT__CK_NORMALIZE_PASSIVE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check the normalization of passive scalars. Make sure to turn on
[[OPT__NORMALIZE_PASSIVE | Runtime-Parameters:-Hydro#OPT__NORMALIZE_PASSIVE ]].
Otherwise this check will likely fail.
    * **Restriction:**

<a name="OPT__CK_RESTRICT"></a>
* #### `OPT__CK_RESTRICT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check the data restriction (i.e., the averages of fine-grid data equal the
coarse-grid data). Make sure to turn on both
[[OPT__FIXUP_RESTRICT | Runtime-Parameters:-Hydro#OPT__FIXUP_RESTRICT]] and
[[OPT__INIT_RESTRICT | Runtime-Parameters:-Initial-Conditions#OPT__INIT_RESTRICT]].
Otherwise this check will likely fail.
    * **Restriction:**

<a name="OPT__CK_FINITE"></a>
* #### `OPT__CK_FINITE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check if all variables are finite (i.e., not "NaN" or "Infinite").
    * **Restriction:**

<a name="OPT__CK_PATCH_ALLOCATE"></a>
* #### `OPT__CK_PATCH_ALLOCATE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check if all patches are properly allocated. This is used mainly for debugging purposes.
    * **Restriction:**

<a name="OPT__CK_FLUX_ALLOCATE"></a>
* #### `OPT__CK_FLUX_ALLOCATE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check if all flux arrays are properly allocated. This is used mainly for debugging purposes.
    * **Restriction:**

<a name="OPT__CK_NEGATIVE"></a>
* #### `OPT__CK_NEGATIVE` &ensp; (0=off, 1=mass density, 2=pressure and entropy, 3=both) &ensp; [0]
    * **Description:**
Check if any of the target variables becomes negative.
    * **Restriction:**

<a name="OPT__CK_MEMFREE"></a>
* #### `OPT__CK_MEMFREE` &ensp; (0=off, >0.0 &#8594; threshold in GB) &ensp; [1.0]
    * **Description:**
Check the free memory. The program will be terminated automatically if the
available memory becomes smaller than a given threshold.
    * **Restriction:**

<a name="OPT__CK_PARTICLE"></a>
* #### `OPT__CK_PARTICLE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check the particle allocation. This is used mainly for debugging purposes.
    * **Restriction:**

<a name="OPT__CK_INTERFACE_B"></a>
* #### `OPT__CK_INTERFACE_B` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check the consistency of the magnetic field on the patch interfaces.
This is used mainly for debugging purposes.
    * **Restriction:**
For [[--mhd | Installation:-Option-List#--mhd]] only.

<a name="OPT__CK_DIVERGENCE_B"></a>
* #### `OPT__CK_DIVERGENCE_B` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check the divergence-free constraint on the magnetic field.
This is used mainly for debugging purposes.
    * **Restriction:**
For [[--mhd | Installation:-Option-List#--mhd]] only.

<a name="OPT__CK_INPUT_FLUID"></a>
* #### `OPT__CK_INPUT_FLUID` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Check the input data of the fluid solver for debugging purposes.
    * **Restriction:**


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
