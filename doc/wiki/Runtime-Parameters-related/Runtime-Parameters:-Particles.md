Parameters described on this page:
[PAR_INIT](#PAR_INIT), &nbsp;
[PAR_NPAR](#PAR_NPAR), &nbsp;
[PAR_IC_FORMAT](#PAR_IC_FORMAT), &nbsp;
[PAR_IC_FLOAT8](#PAR_IC_FLOAT8), &nbsp;
[PAR_IC_INT8](#PAR_IC_INT8), &nbsp;
[PAR_IC_MASS](#PAR_IC_MASS), &nbsp;
[PAR_IC_TYPE](#PAR_IC_TYPE), &nbsp;
[PAR_INTERP](#PAR_INTERP), &nbsp;
[PAR_INTEG](#PAR_INTEG), &nbsp;
[PAR_TR_INTERP](#PAR_TR_INTERP), &nbsp;
[PAR_TR_INTEG](#PAR_TR_INTEG), &nbsp;
[PAR_IMPROVE_ACC](#PAR_IMPROVE_ACC), &nbsp;
[PAR_PREDICT_POS](#PAR_PREDICT_POS), &nbsp;
[PAR_REMOVE_CELL](#PAR_REMOVE_CELL), &nbsp;
[OPT__FREEZE_PAR](#OPT__FREEZE_PAR) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="PAR_INIT"></a>
* #### `PAR_INIT` &ensp; (1=function, 2=restart, 3=`PAR_IC` file) &ensp; [none]
    * **Description:**
Initialization methods for particles.
`PAR_INIT=1`: using a particle initialization function; see
[[Setting IC from Analytical Functions &#8212; Particles | Initial Conditions#IC-Func-Particles]].
`PAR_INIT=2`:
[[restarting from a simulation snapshot | Running the Code#restarting-from-a-snapshot]].
`PAR_INIT=3`: loading a binary file named `PAR_IC`; see
[[Setting IC from Files &#8212; Particles | Initial Conditions#IC-File-Particles]].
    * **Restriction:**
It will be automatically reset to `PAR_INIT=2` when adopting
[[OPT__INIT | Runtime-Parameters:-Initial-Conditions#OPT__INIT]]=2.

<a name="PAR_NPAR"></a>
* #### `PAR_NPAR` &ensp; (&#8805;0) &ensp; [none]
    * **Description:**
Total number of particles.
    * **Restriction:**
Must be set when [PAR_INIT](#PAR_INIT)=1/3;
must be consistent with the number of particles loaded from an initial
condition file (e.g., for [PAR_INIT](#PAR_INIT)=3 and test problems
with their own particle initializers);
useless when restarting from a snapshot (i.e., [PAR_INIT](#PAR_INIT)=2);
**must be an integer**.

<a name="PAR_IC_FORMAT"></a>
* #### `PAR_IC_FORMAT` &ensp; (1=[attribute][id], 2=[id][attribute]; row-major) &ensp; [1]
    * **Description:**
Data format of the particle initial condition file `PAR_IC` used by
[PAR_INIT](#PAR_INIT)=3. See
[[Setting IC from Files &#8212; Particles | Initial Conditions#IC-File-Particles]]
for details.
    * **Restriction:**

<a name="PAR_IC_FLOAT8"></a>
* #### `PAR_IC_FLOAT8` &ensp; (<0: same as [[--double_par | Installation:-Option-List#--double_par]], 0=single precision, 1=double precision) &ensp; [-1]
    * **Description:**
Floating-point precision of the particle initial condition file `PAR_IC`.
    * **Restriction:**

<a name="PAR_IC_INT8"></a>
* #### `PAR_IC_INT8` &ensp; (<0: same as [[--long_par | Installation:-Option-List#--long_par]], 0=32-bit integer (`int`), 1=64-bit integer (`long`)) &ensp; [-1]
    * **Description:**
Integer width of the particle initial condition file `PAR_IC`.
    * **Restriction:**

<a name="PAR_IC_MASS"></a>
* #### `PAR_IC_MASS` &ensp; (&#8805;0.0; <0.0 &#8594; off) &ensp; [-1.0]
    * **Description:**
Assigning this particle mass to all particles when adopting
[PAR_INIT](#PAR_INIT)=3. Note that when enabling this functionality (by setting `PAR_IC_MASS`&#8805;0.0),
the particle initial condition file `PAR_IC` should not include the particle mass data.
See also
[[Setting IC from Files &#8212; Particles | Initial Conditions#IC-File-Particles]].
    * **Restriction:**

<a name="PAR_IC_TYPE"></a>
* #### `PAR_IC_TYPE` &ensp; (&#8805;0; <0 &#8594; off) &ensp; [-1]
    * **Description:**
Assigning this particle type to all particles when adopting
[PAR_INIT](#PAR_INIT)=3. Note that when enabling this functionality (by setting `PAR_IC_TYPE`&#8805;0),
the particle initial condition file `PAR_IC` should not include the particle type data.
See also
[[Setting IC from Files &#8212; Particles | Initial Conditions#IC-File-Particles]].
    * **Restriction:**

<a name="PAR_INTERP"></a>
* #### `PAR_INTERP` &ensp; (1=NGP, 2=CIC, 3=TSC) &ensp; [2]
    * **Description:**
Massive particle interpolation schemes:
NGP=nearest-grid-point, CIC=cloud-in-cell, TSC=triangular-shape cloud.
Accuracy: TSC > CIC > NGP. Performance: NGP > CIC > TSC. In general,
NGP is not recommended.
    * **Restriction:**

<a name="PAR_INTEG"></a>
* #### `PAR_INTEG` &ensp; (1=Euler, 2=KDK) &ensp; [2]
    * **Description:**
Massive particle integration scheme. Euler integration is only first-order accurate
and is generally not recommended.
    * **Restriction:**

<a name="PAR_TR_INTERP"></a>
* #### `PAR_TR_INTERP` &ensp; (1=NGP, 2=CIC, 3=TSC) &ensp; [3]
    * **Description:**
Tracer particle interpolation schemes. See [PAR_INTERP](#PAR_INTERP) for details.
    * **Restriction:**

<a name="PAR_TR_INTEG"></a>
* #### `PAR_TR_INTEG` &ensp; (1=Euler, 2=RK2) &ensp; [2]
    * **Description:**
Tracer particle integration scheme. Euler integration is only first-order accurate
and is generally not recommended.
    * **Restriction:**

<a name="PAR_IMPROVE_ACC"></a>
* #### `PAR_IMPROVE_ACC` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Improve the force accuracy for particles close to the patch boundaries.
    * **Restriction:**
Only applicable when adopting [PAR_INTERP](#PAR_INTERP)=2/3 and
enabling the compilation option [[--store_pot_ghost | Installation:-Option-List#--store_pot_ghost]].

<a name="PAR_PREDICT_POS"></a>
* #### `PAR_PREDICT_POS` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Predict particle position during mass assignment to
synchronize particles and grids.
    * **Restriction:**

<a name="PAR_REMOVE_CELL"></a>
* #### `PAR_REMOVE_CELL` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [depend]
    * **Description:**
Remove particles too close to the simulation boundaries when adopting
non-periodic boundary conditions. It is necessary because the force
accuracy around the simulation boundaries may be deteriorated due
to the potential extrapolation. Specifically, we remove particles
located in the outermost `PAR_REMOVE_CELL` root cells in the simulation
box. The default value is 1.0 ~ 2.0 depending on the adopted particle
interpolation scheme ([PAR_INTERP](#PAR_INTERP)).
    * **Restriction:**
Only applicable when adopting the isolated gravity boundary condition
(i.e., [[OPT__BC_POT | Runtime-Parameters:-Gravity#OPT__BC_POT]]=2).

<a name="OPT__FREEZE_PAR"></a>
* #### `OPT__FREEZE_PAR` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Do not update particle position and velocity (except for tracer particles).
It can be useful for evolving fluid in a static gravitational potential of particles.
    * **Restriction:**


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Particles | Particles]]
