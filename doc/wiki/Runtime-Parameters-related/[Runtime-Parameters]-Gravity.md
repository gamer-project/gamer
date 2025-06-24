Parameters described on this page:
[OPT__BC_POT](#OPT__BC_POT), &nbsp;
[GFUNC_COEFF0](#GFUNC_COEFF0), &nbsp;
[NEWTON_G](#NEWTON_G), &nbsp;
[SOR_OMEGA](#SOR_OMEGA), &nbsp;
[SOR_MAX_ITER](#SOR_MAX_ITER), &nbsp;
[SOR_MIN_ITER](#SOR_MIN_ITER), &nbsp;
[MG_MAX_ITER](#MG_MAX_ITER), &nbsp;
[MG_NPRE_SMOOTH](#MG_NPRE_SMOOTH), &nbsp;
[MG_NPOST_SMOOTH](#MG_NPOST_SMOOTH), &nbsp;
[MG_TOLERATED_ERROR](#MG_TOLERATED_ERROR), &nbsp;
[OPT__GRA_P5_GRADIENT](#OPT__GRA_P5_GRADIENT), &nbsp;
[OPT__SELF_GRAVITY](#OPT__SELF_GRAVITY), &nbsp;
[OPT__EXT_ACC](#OPT__EXT_ACC), &nbsp;
[OPT__EXT_POT](#OPT__EXT_POT), &nbsp;
[EXT_POT_TABLE_NAME](#EXT_POT_TABLE_NAME), &nbsp;
[EXT_POT_TABLE_NPOINT_X](#EXT_POT_TABLE_NPOINT_X), &nbsp;
[EXT_POT_TABLE_NPOINT_Y](#EXT_POT_TABLE_NPOINT_Y), &nbsp;
[EXT_POT_TABLE_NPOINT_Z](#EXT_POT_TABLE_NPOINT_Z), &nbsp;
[EXT_POT_TABLE_DH](#EXT_POT_TABLE_DH), &nbsp;
[EXT_POT_TABLE_EDGEL_X](#EXT_POT_TABLE_EDGEL_X), &nbsp;
[EXT_POT_TABLE_EDGEL_Y](#EXT_POT_TABLE_EDGEL_Y), &nbsp;
[EXT_POT_TABLE_EDGEL_Z](#EXT_POT_TABLE_EDGEL_Z) &nbsp;



Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

 <a name="OPT__BC_POT"></a>
* #### `OPT__BC_POT` &ensp; (1=periodic, 2=isolated) &ensp; [none]
    * **Description:**
Gravity boundary condition. See also
[[Potential Outside the Isolated Boundaries | Runtime-Parameters:-Refinement#potential-outside-the-isolated-boundaries]].
    * **Restriction:**

<a name="GFUNC_COEFF0"></a>
* #### `GFUNC_COEFF0` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [depend]
    * **Description:**
Green's function coefficient for calculating the "self"
gravitational potential of each cell (i.e., the gravitational
potential contribution from the mass within the same cell).
The default value depends on the particle interpolation scheme
(`PAR_INTERP`). Set GFUNC_COEFF0=0.0 to exclude the self-potential.
    * **Restriction:**
Only applicable to the isolated boundary condition.

<a name="NEWTON_G"></a>
* #### `NEWTON_G` &ensp; (>0.0) &ensp; [conform to the unit system set by [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]]]
    * **Description:**
Gravitational constant.
    * **Restriction:**
It will be overwritten by the default value when [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]]
is on; no default when [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]] is off.

<a name="SOR_OMEGA"></a>
* #### `SOR_OMEGA` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [1.69]
    * **Description:**
Parameter of the SOR Poisson solver for optimizing the convergence rate.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=SOR.

<a name="SOR_MAX_ITER"></a>
* #### `SOR_MAX_ITER` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [single precision=60, double precision=110]
    * **Description:**
Maximum number of iterations in the SOR Poisson solver.
The default value depends on the adopted floating-point accuracy (`FLOAT8`).
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=SOR.

<a name="SOR_MIN_ITER"></a>
* #### `SOR_MIN_ITER` &ensp; (&#8805;3; <0 &#8594; set to default) &ensp; [10]
    * **Description:**
Minimum number of iterations in the SOR Poisson solver.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=SOR.

<a name="MG_MAX_ITER"></a>
* #### `MG_MAX_ITER` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [single precision=10, double precision=20]
    * **Description:**
Maximum number of iterations in the multigrid Poisson solver.
The default value depends on the adopted floating-point accuracy ([[--double | Installation:-Option-List#--double]]).
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=MG.

<a name="MG_NPRE_SMOOTH"></a>
* #### `MG_NPRE_SMOOTH` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [3]
    * **Description:**
Number of pre-smoothing steps in the multigrid Poisson solver.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=MG.

<a name="MG_NPOST_SMOOTH"></a>
* #### `MG_NPOST_SMOOTH` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [3]
    * **Description:**
Number of post-smoothing steps in the multigrid Poisson solver.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=MG.

<a name="MG_TOLERATED_ERROR"></a>
* #### `MG_TOLERATED_ERROR` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [single precision=1e-6, double precision=1e-15]
    * **Description:**
Maximum tolerable error in the multigrid Poisson solver.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=MG.

<a name="OPT__GRA_P5_GRADIENT"></a>
* #### `OPT__GRA_P5_GRADIENT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Use 5-point instead of 3-point gradient for calculating the
gravitational acceleration from potential.
*This functionality is only experimental.*
    * **Restriction:**
Must manually set `#define GRA_GHOST_SIZE 2` (and `#define USG_GHOST_SIZE 2`
as well when adopting the compilation option
[[--unsplit_gravity | Installation:-Option-List#--unsplit_gravity]])
in the header file `Macro.h`. Unsupported for particle update.

<a name="OPT__SELF_GRAVITY"></a>
* #### `OPT__SELF_GRAVITY` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Enable self-gravity.
    * **Restriction:**

<a name="OPT__EXT_ACC"></a>
* #### `OPT__EXT_ACC` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Add external acceleration. See
[External Acceleration/Potential](#external-accelerationpotential)
for how to specify external acceleration.
    * **Restriction:**
Not applicable to the wave dark matter simulations
([[--model | Installation:-Option-List#--model]]=ELBDM).

<a name="OPT__EXT_POT"></a>
* #### `OPT__EXT_POT` &ensp; (0=off, 1=function, 2=table) &ensp; [0]
    * **Description:**
Add external potential. See
[External Acceleration/Potential](#external-accelerationpotential)
for how to specify external potential.
    * **Restriction:**

<a name="EXT_POT_TABLE_NAME"></a>
* #### `EXT_POT_TABLE_NAME` &ensp; (string) &ensp; [none]
    * **Description:**
For [OPT__EXT_POT](#OPT__EXT_POT)`=2`: table name.
    * **Restriction:**

<a name="EXT_POT_TABLE_NPOINT_X"></a>
* #### `EXT_POT_TABLE_NPOINT_X` &ensp; (>1) &ensp; [none]
    * **Description:**
For [OPT__EXT_POT](#OPT__EXT_POT)`=2`: table size (i.e., number of data points) along the x direction.
    * **Restriction:**

<a name="EXT_POT_TABLE_NPOINT_Y"></a>
* #### `EXT_POT_TABLE_NPOINT_Y` &ensp; (>1) &ensp; [none]
    * **Description:**
See [EXT_POT_TABLE_NPOINT_X](#EXT_POT_TABLE_NPOINT_X).
    * **Restriction:**

<a name="EXT_POT_TABLE_NPOINT_Z"></a>
* #### `EXT_POT_TABLE_NPOINT_Z` &ensp; (>1) &ensp; [none]
    * **Description:**
See [EXT_POT_TABLE_NPOINT_X](#EXT_POT_TABLE_NPOINT_X).
    * **Restriction:**

<a name="EXT_POT_TABLE_DH"></a>
* #### `EXT_POT_TABLE_DH` &ensp; (>0.0) &ensp; [none]
    * **Description:**
For [OPT__EXT_POT](#OPT__EXT_POT)`=2`: spatial interval between adjacent tabular data.
    * **Restriction:**

<a name="EXT_POT_TABLE_EDGEL_X"></a>
* #### `EXT_POT_TABLE_EDGEL_X` &ensp; (cover the simulation domain plus ghost zones) &ensp; [none]
    * **Description:**
For [OPT__EXT_POT](#OPT__EXT_POT)`=2`: starting x coordinate of the tabular data.
    * **Restriction:**
Must cover the simulation domain plus gravity ghost zones. Specifically, `EXT_POT_TABLE_EDGEL_X` must be
smaller or equal to `DomainLeftEdge - 1*RootLevelCell` and
`EXT_POT_TABLE_EDGEL_X + (EXT_POT_TABLE_NPOINT_X-1)*EXT_POT_TABLE_DH`
must be larger or equal to `DomainRightEdge + 1*RootLevelCell`.

<a name="EXT_POT_TABLE_EDGEL_Y"></a>
* #### `EXT_POT_TABLE_EDGEL_Y` &ensp; (cover the simulation domain) &ensp; [none]
    * **Description:**
See [EXT_POT_TABLE_EDGEL_X](#EXT_POT_TABLE_EDGEL_X).
    * **Restriction:**

<a name="EXT_POT_TABLE_EDGEL_Z"></a>
* #### `EXT_POT_TABLE_EDGEL_Z` &ensp; (cover the simulation domain) &ensp; [none]
    * **Description:**
See [EXT_POT_TABLE_EDGEL_X](#EXT_POT_TABLE_EDGEL_X).
    * **Restriction:**


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Gravity | Gravity]]
