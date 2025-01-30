Parameters described on this page:
[OPT__INIT](#OPT__INIT), &nbsp;
[OPT__INIT_BFIELD_BYVECPOT](#OPT__INIT_BFIELD_BYVECPOT), &nbsp;
[RESTART_LOAD_NRANK](#RESTART_LOAD_NRANK), &nbsp;
[OPT__RESTART_RESET](#OPT__RESTART_RESET), &nbsp;
[OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL), &nbsp;
[OPT__UM_IC_NLEVEL](#OPT__UM_IC_NLEVEL), &nbsp;
[OPT__UM_IC_NVAR](#OPT__UM_IC_NVAR), &nbsp;
[OPT__UM_IC_FORMAT](#OPT__UM_IC_FORMAT), &nbsp;
[OPT__UM_IC_FLOAT8](#OPT__UM_IC_FLOAT8), &nbsp;
[OPT__UM_IC_DOWNGRADE](#OPT__UM_IC_DOWNGRADE), &nbsp;
[OPT__UM_IC_REFINE](#OPT__UM_IC_REFINE), &nbsp;
[OPT__UM_IC_LOAD_NRANK](#OPT__UM_IC_LOAD_NRANK), &nbsp;
[OPT__INIT_RESTRICT](#OPT__INIT_RESTRICT), &nbsp;
[INIT_SUBSAMPLING_NCELL](#INIT_SUBSAMPLING_NCELL), &nbsp;
[OPT__FFTW_STARTUP](#OPT__FFTW_STARTUP) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="OPT__INIT"></a>
* #### `OPT__INIT` &ensp; (1=function, 2=restart, 3=`UM_IC` file) &ensp; [none]
    * **Description:**
Grid initialization method.
`OPT__INIT=1`: using analytical functions; see
[[Setting IC from Analytical Functions &#8212; Grids | Initial-Conditions#IC-Func-Grids]].
`OPT__INIT=2`:
[[restarting from a simulation snapshot | Running-the-Code#restarting-from-a-snapshot]].
`OPT__INIT=3`: loading a uniform-mesh binary file named `UM_IC`; see
[[Setting IC from Files &#8212; Grids | Initial-Conditions#IC-File-Grids]].
    * **Restriction:**

<a name="OPT__INIT_BFIELD_BYVECPOT"></a>
* #### `OPT__INIT_BFIELD_BYVECPOT` &ensp; (0=off, 1=file, 2=function) &ensp; [0]
    * **Description:**
Set the magnetic field from either a vector potential file named `B_IC` (see
[[Setting IC from Files &#8212; Magnetic Field | Initial-Conditions#IC-File-BField]]) or an
analytical vector potential function (see
[[Setting IC from Functions &#8212; Magnetic Field | Initial-Conditions#IC-Func-BField]]).
    * **Restriction:**
For [[--mhd | Installation:-Option-List#--mhd]] only.

<a name="RESTART_LOAD_NRANK"></a>
* #### `RESTART_LOAD_NRANK` &ensp; (>0) &ensp; [1]
    * **Description:**
Number of parallel I/O for restart. In other words, `RESTART_LOAD_NRANK`
MPI processes will load the restart file in parallel.
    * **Restriction:**

<a name="OPT__RESTART_RESET"></a>
* #### `OPT__RESTART_RESET` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
During restart, reset some of the simulation status parameters
(e.g., step, time, snapshot ID) to their initial values as if the
simulation starts over again.
    * **Restriction:**

<a name="OPT__UM_IC_LEVEL"></a>
* #### `OPT__UM_IC_LEVEL` &ensp; (0 &#8804; input < [[--nlevel | Installation:-Option-List#--nlevel]]) &ensp; [0]
    * **Description:**
Starting AMR level in the uniform-mesh initial condition file.
See [[Setting IC from Files &#8212; Grids | Initial-Conditions#IC-File-Grids]] for details.
    * **Restriction:**

<a name="OPT__UM_IC_NLEVEL"></a>
* #### `OPT__UM_IC_NLEVEL` &ensp; (1 &#8804; input &#8804; [[--nlevel | Installation:-Option-List#--nlevel]]-[OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL)) &ensp; [1]
    * **Description:**
Number of AMR levels in the uniform-mesh initial condition file.
See [[Setting IC from Files &#8212; Grids | Initial-Conditions#IC-File-Grids]] for details.
    * **Restriction:**

<a name="OPT__UM_IC_NVAR"></a>
* #### `OPT__UM_IC_NVAR` &ensp; (0 &#8804; input < total number of gas fields) &ensp; [total number of gas fields]
    * **Description:**
Number of fluid variables stored in the uniform-mesh initial condition file.
The default value is
5+[[--passive | Installation:-Option-List#--passive]]
for [[--model | Installation:-Option-List#--model]]=HYDRO.
See [[Setting IC from Files &#8212; Grids | Initial-Conditions#IC-File-Grids]] for details.
    * **Restriction:**

<a name="OPT__UM_IC_FORMAT"></a>
* #### `OPT__UM_IC_FORMAT` &ensp; (1=[v][z][y][x], 2=[z][y][x][v]; row-major and v=field) &ensp; [1]
    * **Description:**
Data format of the uniform-mesh initial condition file. See
[[Setting IC from Files &#8212; Grids | Initial-Conditions#IC-File-Grids]] for details.
    * **Restriction:**

<a name="OPT__UM_IC_FLOAT8"></a>
* #### `OPT__UM_IC_FLOAT8` &ensp; (<0: same as [[--double | Installation:-Option-List#--double]], 0=single precision, 1=double precision) &ensp; [-1]
    * **Description:**
Floating-point precision of the uniform-mesh initial condition file.
    * **Restriction:**

<a name="OPT__UM_IC_DOWNGRADE"></a>
* #### `OPT__UM_IC_DOWNGRADE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Downgrade (i.e. derefine) the uniform-mesh initial condition data for cells
not satisfying any refinement criteria.
See [[Setting IC from Files &#8212; Grids | Initial-Conditions#IC-File-Grids]] for details.
    * **Restriction:**

<a name="OPT__UM_IC_REFINE"></a>
* #### `OPT__UM_IC_REFINE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Refine the uniform-mesh initial condition data from level `OPT__UM_IC_LEVEL` to
[[MAX_LEVEL | Runtime Parameters:-Refinement#MAX_LEVEL]] for cells satisfying the adopted
refinement criteria.
See [[Setting IC from Files &#8212; Grids | Initial-Conditions#IC-File-Grids]] for details.
    * **Restriction:**

<a name="OPT__UM_IC_LOAD_NRANK"></a>
* #### `OPT__UM_IC_LOAD_NRANK` &ensp; (>0) &ensp; [1]
    * **Description:**
Number of parallel I/O for loading the uniform-mesh initial condition file.
Specifically, it allows `OPT__UM_IC_LOAD_NRANK` MPI processes to load the
initial condition file concurrently. But the actually achieved parallel I/O
depends on the system specifications.
See also [[Setting IC from Files &#8212; Grids | Initial-Conditions#IC-File-Grids]].
    * **Restriction:**

<a name="OPT__INIT_RESTRICT"></a>
* #### `OPT__INIT_RESTRICT` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
For non-leaf patches, replace fluid data by the volume-weighted average
of their child patch data. It is similar to the option
[[OPT__FIXUP_RESTRICT | Runtime-Parameters:-Hydro#OPT__FIXUP_RESTRICT]]
except that it only applies to the initial condition.
    * **Restriction:**

<a name="INIT_SUBSAMPLING_NCELL"></a>
* #### `INIT_SUBSAMPLING_NCELL` &ensp; (0=off, >0 &#8594; number of sub-cells along each direction) &ensp; [0]
    * **Description:**
Perform sub-sampling when constructing the grid IC to make it smoother.
Specifically, each cell will be divided into <var>N</var><sub>sub</sub><sup>3</sup>
sub-cells when calling the grid IC function, where
<var>N</var><sub>sub</sub> = `INIT_SUBSAMPLING_NCELL`, and then
take the volume-weighted average of these sub-cells.
    * **Restriction:**
Only applicable when adopting [OPT__INIT](#OPT__INIT)=1.

<a name="OPT__FFTW_STARTUP"></a>
* #### `OPT__FFTW_STARTUP` &ensp; (-1 &#8594; set to default, 0=ESTIMATE, 1=MEASURE, 2=PATIENT) &ensp; [-1]
    * **Description:**
Initialize FFTW plans. `MEASURE` is recommended for the balance
between FFTW plan initialization time and FFT performance.
Note that simulation results can vary in each run on the level of
machine precision for `OPT__FFTW_STARTUP != ESTIMATE`.
    * **Restriction:**
`PATIENT` is not supported by FFTW2.
Must use `ESTIMATE` when enabling
[[--bitwise_reproducibility | Installation:-Option-List#--bitwise_reproducibility]].


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Initial Conditions | Initial-Conditions]]
