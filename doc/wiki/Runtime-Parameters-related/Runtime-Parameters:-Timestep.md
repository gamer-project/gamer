This page describes various timestep constraints. See also Sections 2.1 and 2.6 in the
[GAMER-2 code paper](https://arxiv.org/abs/1712.07070).


## Compilation Options

Related options: none


## Runtime Parameters

Parameters described on this page:
[DT__MAX](#DT__MAX), &nbsp;
[DT__FLUID](#DT__FLUID), &nbsp;
[DT__FLUID_INIT](#DT__FLUID_INIT), &nbsp;
[DT__SPEED_OF_LIGHT](#DT__SPEED_OF_LIGHT), &nbsp;
[DT__GRAVITY](#DT__GRAVITY), &nbsp;
[DT__PARVEL](#DT__PARVEL), &nbsp;
[DT__PARVEL_MAX](#DT__PARVEL_MAX), &nbsp;
[DT__PARACC](#DT__PARACC), &nbsp;
[DT__CR_DIFFUSION](#DT__CR_DIFFUSION), &nbsp;
[DT__MAX_DELTA_A](#DT__MAX_DELTA_A), &nbsp;
[DT__SYNC_PARENT_LV](#DT__SYNC_PARENT_LV), &nbsp;
[DT__SYNC_CHILDREN_LV](#DT__SYNC_CHILDREN_LV), &nbsp;
[OPT__DT_USER](#OPT__DT_USER), &nbsp;
[OPT__DT_LEVEL](#OPT__DT_LEVEL), &nbsp;
[OPT__RECORD_DT](#OPT__RECORD_DT), &nbsp;
[AUTO_REDUCE_DT](#AUTO_REDUCE_DT), &nbsp;
[AUTO_REDUCE_DT_FACTOR](#AUTO_REDUCE_DT_FACTOR), &nbsp;
[AUTO_REDUCE_DT_FACTOR_MIN](#AUTO_REDUCE_DT_FACTOR_MIN), &nbsp;
[AUTO_REDUCE_MINMOD_FACTOR](#AUTO_REDUCE_MINMOD_FACTOR), &nbsp;
[AUTO_REDUCE_MINMOD_MIN](#AUTO_REDUCE_MINMOD_MIN), &nbsp;
[AUTO_REDUCE_INT_MONO_FACTOR](#AUTO_REDUCE_INT_MONO_FACTOR), &nbsp;
[AUTO_REDUCE_INT_MONO_MIN](#AUTO_REDUCE_INT_MONO_MIN) &nbsp;


Other related parameters:
[[]] &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="DT__MAX"></a>
* #### `DT__MAX` &ensp; (&#8805;0.0; <0.0 &#8594; off) &ensp; [-1.0]
    * **Description:**
Maximum allowed time-step on all levels. This can be used to, for example, avoid
too large time-steps on lower levels.
    * **Restriction:**

<a name="DT__FLUID"></a>
* #### `DT__FLUID` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [-1.0]
    * **Description:**
Courant–Friedrichs–Lewy (CFL) safety factor for the hydrodynamic solver.
The default value and stable regime depend on the adopted
[[fluid scheme | Installation:-Option-List#--flu_scheme]].
See Section 2.6, Eqs. [1-2] in the [GAMER-2 code paper](https://arxiv.org/abs/1712.07070)
for the exact formulae.
    * **Restriction:**

<a name="DT__FLUID_INIT"></a>
* #### `DT__FLUID_INIT` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [[DT__FLUID](#DT__FLUID)]
    * **Description:**
CFL safety factor for the hydrodynamic solver _at the first step_. This could be
useful when the first step requires a much smaller timestep.
    * **Restriction:**
Useless for restart.

<a name="DT__SPEED_OF_LIGHT"></a>
* #### `DT__SPEED_OF_LIGHT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Using the speed of light as the maximum information propagation speed for determining timesteps
in special relativistic hydrodynamic simulations.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--srhd | Installation:-Option-List#--srhd]].

<a name="DT__GRAVITY"></a>
* #### `DT__GRAVITY` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [-1.0]
    * **Description:**
Safety factor when determining timestep from the gravitational acceleration of fluid.
See Section 2.6, Eq. [3] in the [GAMER-2 code paper](https://arxiv.org/abs/1712.07070)
for the exact formula.
    * **Restriction:**

<a name="DT__PARVEL"></a>
* #### `DT__PARVEL` &ensp; (&#8805;0.0) &ensp; [0.5]
    * **Description:**
Safety factor when determining timestep from the particle velocity.
See Section 2.6, Eq. [4] in the [GAMER-2 code paper](https://arxiv.org/abs/1712.07070)
for the exact formula.
    * **Restriction:**

<a name="DT__PARVEL_MAX"></a>
* #### `DT__PARVEL_MAX` &ensp; (&#8805;0.0; <0.0 &#8594; off) &ensp; [-1.0]
    * **Description:**
Maximum allowed value of the timestep determined from the particle velocity.
This could be useful when, for example, all particles are initially at rest,
for which the timestep determined from the particle velocity is infinity.
    * **Restriction:**

<a name="DT__PARACC"></a>
* #### `DT__PARACC` &ensp; (>0.0; &#8804;0.0 &#8594; off) &ensp; [0.5]
    * **Description:**
Safety factor when determining timestep from the particle acceleration.
See Section 2.6, Eq. [3] in the [GAMER-2 code paper](https://arxiv.org/abs/1712.07070)
for the exact formula.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--store_par_acc | Installation:-Option-List#--store_par_acc]].

<a name="DT__CR_DIFFUSION"></a>
* #### `DT__CR_DIFFUSION` &ensp; (&#8805;0.0) &ensp; [0.3]
    * **Description:**
CFL safety factor for cosmic-ray diffusion.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--cr_diffusion | Installation:-Option-List#--cr_diffusion]].

<a name="DT__MAX_DELTA_A"></a>
* #### `DT__MAX_DELTA_A` &ensp; (&#8805;0.0) &ensp; [0.01]
    * **Description:**
Maximum allowed fraction of increase in the cosmic scale factor <var>a</var>.
Specifically, it enforces &Delta;<var>a</var> &#8804; `DT__MAX_DELTA_A` &#8901; <var>a</var>.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--comoving | Installation:-Option-List#--comoving]].

<a name="DT__SYNC_PARENT_LV"></a>
* #### `DT__SYNC_PARENT_LV` &ensp; (&#8805;0.0) &ensp; [0.1]
    * **Description:**
Allow timestep to _increase_ by `1.0+DT__SYNC_PARENT` to help synchronize
with the parent level. See also Section 2.1 in the
[GAMER-2 code paper](https://arxiv.org/abs/1712.07070) for more details
    * **Restriction:**
For [OPT__DT_LEVEL](#OPT__DT_LEVEL)=3 only.

<a name="DT__SYNC_CHILDREN_LV"></a>
* #### `DT__SYNC_CHILDREN_LV` &ensp; (&#8805;0.0; <0.0 &#8594; off) &ensp; [0.1]
    * **Description:**
Allow timestep to _decrease_ by `1.0+DT__SYNC_CHILDREN_LV` to help synchronize
with the children level. See also Section 2.1 in the
[GAMER-2 code paper](https://arxiv.org/abs/1712.07070) for more details
    * **Restriction:**
For [OPT__DT_LEVEL](#OPT__DT_LEVEL)=3 only.

<a name="OPT__DT_USER"></a>
* #### `OPT__DT_USER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Enable a user-defined timestep criterion. Edit `src/Miscellaneous/Mis_GetTimeStep_User.cpp`
or a problem-specific function (for the latter, see
[[Add Problem Specific Functionalities | Adding-New-Simulations#vi-add-problem-specific-functionalities]]).
    * **Restriction:**

<a name="OPT__DT_LEVEL"></a>
* #### `OPT__DT_LEVEL` &ensp; (1=shared, 2=differ by two, 3=flexible) &ensp; [3]
    * **Description:**
Constraints on the timesteps of adjacent levels.
        * `1`: All levels must share the same timestep
        * `2`: Timesteps of a parent level is fixed to be twice smaller than its parent level
        * `3`: No constraint (except that the timestep of a children level cannot be larger than that of a parent level)
    * **Restriction:**

<a name="OPT__RECORD_DT"></a>
* #### `OPT__RECORD_DT` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Record the timesteps set by various constraints in the file
[[Record__TimeStep | Simulation-Logs:-Record__TimeStep]].
    * **Restriction:**

<a name="AUTO_REDUCE_DT"></a>
* #### `AUTO_REDUCE_DT` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Automatically reduce timestep by a factor of
[AUTO_REDUCE_DT_FACTOR](#AUTO_REDUCE_DT_FACTOR),
[[MINMOD_COEFF | Runtime-Parameters:-Hydro#MINMOD_COEFF]] by a factor of
[AUTO_REDUCE_MINMOD_FACTOR](#AUTO_REDUCE_MINMOD_FACTOR), and
[[INT_MONO_COEFF | Runtime-Parameters:-Interpolation#INT_MONO_COEFF]]
by a factor of
[AUTO_REDUCE_INT_MONO_FACTOR](#AUTO_REDUCE_INT_MONO_FACTOR)
when the program fails until the timestep reducing factor
becomes smaller than
[AUTO_REDUCE_DT_FACTOR_MIN](#AUTO_REDUCE_DT_FACTOR_MIN)
or the interpolation coefficients become smaller than
[AUTO_REDUCE_MINMOD_MIN](#AUTO_REDUCE_MINMOD_MIN) or
[AUTO_REDUCE_INT_MONO_MIN](#AUTO_REDUCE_INT_MONO_MIN).
    * **Restriction:**
For [OPT__DT_LEVEL](#OPT__DT_LEVEL)=3 only.

<a name="AUTO_REDUCE_DT_FACTOR"></a>
* #### `AUTO_REDUCE_DT_FACTOR` &ensp; (>0.0) &ensp; [1.0]
    * **Description:**
See [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).
    * **Restriction:**

<a name="AUTO_REDUCE_DT_FACTOR_MIN"></a>
* #### `AUTO_REDUCE_DT_FACTOR_MIN` &ensp; (&#8805;0.0) &ensp; [0.1]
    * **Description:**
See [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).
    * **Restriction:**

<a name="AUTO_REDUCE_MINMOD_FACTOR"></a>
* #### `AUTO_REDUCE_MINMOD_FACTOR` &ensp; (>0.0) &ensp; [0.8]
    * **Description:**
See [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).
    * **Restriction:**
Must enable [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).

<a name="AUTO_REDUCE_MINMOD_MIN"></a>
* #### `AUTO_REDUCE_MINMOD_MIN` &ensp; (&#8805;0.0) &ensp; [0.01]
    * **Description:**
See [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).
    * **Restriction:**
Must enable [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).

<a name="AUTO_REDUCE_INT_MONO_FACTOR"></a>
* #### `AUTO_REDUCE_INT_MONO_FACTOR` &ensp; (>0.0) &ensp; [0.8]
    * **Description:**
See [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).
    * **Restriction:**
Must enable [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).

<a name="AUTO_REDUCE_INT_MONO_MIN"></a>
* #### `AUTO_REDUCE_INT_MONO_MIN` &ensp; (&#8805;0.0) &ensp; [0.01]
    * **Description:**
See [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).
    * **Restriction:**
Must enable [AUTO_REDUCE_DT](#AUTO_REDUCE_DT).


<br>

## Links
* [[Main page of Runtime Parameters | Runtime-Parameters]]
