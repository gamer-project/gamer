This page describes the supported interpolation schemes for (i) grid refinement
and (ii) preparing the ghost zones around any coarse-fine boundaries.


## Compilation Options

Related options: None


## Runtime Parameters

Parameters described on this page:
[OPT__INT_TIME](#OPT__INT_TIME), &nbsp;
[OPT__INT_PRIM](#OPT__INT_PRIM), &nbsp;
[OPT__FLU_INT_SCHEME](#OPT__FLU_INT_SCHEME), &nbsp;
[OPT__REF_FLU_INT_SCHEME](#OPT__REF_FLU_INT_SCHEME), &nbsp;
[OPT__MAG_INT_SCHEME](#OPT__MAG_INT_SCHEME), &nbsp;
[OPT__REF_MAG_INT_SCHEME](#OPT__REF_MAG_INT_SCHEME), &nbsp;
[OPT__POT_INT_SCHEME](#OPT__POT_INT_SCHEME), &nbsp;
[OPT__RHO_INT_SCHEME](#OPT__RHO_INT_SCHEME), &nbsp;
[OPT__GRA_INT_SCHEME](#OPT__GRA_INT_SCHEME), &nbsp;
[OPT__REF_POT_INT_SCHEME](#OPT__REF_POT_INT_SCHEME), &nbsp;
[INT_MONO_COEFF](#INT_MONO_COEFF), &nbsp;
[INT_MONO_COEFF_B](#INT_MONO_COEFF_B), &nbsp;
[MONO_MAX_ITER](#MONO_MAX_ITER), &nbsp;
[INT_OPP_SIGN_0TH_ORDER](#INT_OPP_SIGN_0TH_ORDER) &nbsp;

Other related parameters:
[[AUTO_REDUCE_INT_MONO_FACTOR | Runtime-Parameters:-Timestep#AUTO_REDUCE_INT_MONO_FACTOR]] &nbsp;

<a name="INT_TABLE"></a>
Supported interpolation schemes:

| ID | Interpolation Scheme|
|---:|:---|
|-1 | Set to default|
|1 | 3D MinMod limiter|
|2 | 1D MinMod limiter|
|3 | vanLeer limiter|
|4 | Conservative quadratic|
|5 | Non-conservative quadratic|
|6 | Conservative quartic|
|7 | Non-conservative quartic|


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**


<a name="OPT__INT_TIME"></a>
* #### `OPT__INT_TIME` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Enable temporal interpolation when preparing ghost zones.
    * **Restriction:**
Only applicable when adopting the adaptive timestep integration
(i.e., [[OPT__DT_LEVEL | Runtime-Parameters:-Timestep#OPT__DT_LEVEL]]=2/3).

<a name="OPT__INT_PRIM"></a>
* #### `OPT__INT_PRIM` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Perform interpolation on primitive variables when interpolation
on conserved variables fails.
    * **Restriction:**

<a name="OPT__FLU_INT_SCHEME"></a>
* #### `OPT__FLU_INT_SCHEME` &ensp; (see ["Supported interpolation schemes"](#INT_TABLE)) &ensp; [-1]
    * **Description:**
Interpolation scheme for preparing the ghost zones of the fluid solver.
    * **Restriction:**

<a name="OPT__REF_FLU_INT_SCHEME"></a>
* #### `OPT__REF_FLU_INT_SCHEME` &ensp; (see ["Supported interpolation schemes"](#INT_TABLE)) &ensp; [-1]
    * **Description:**
Interpolation scheme for computing the fluid variables on the newly refined patches.
    * **Restriction:**

<a name="OPT__MAG_INT_SCHEME"></a>
* #### `OPT__MAG_INT_SCHEME` &ensp; (only 2, 3, 4, 6 in ["Supported interpolation schemes"](#INT_TABLE)) &ensp; [4]
    * **Description:**
Interpolation scheme for preparing the ghost-zone magnetic field of the fluid solver.
    * **Restriction:**
For [[--mhd | Installation:-Option-List#--mhd]] only.

<a name="OPT__REF_MAG_INT_SCHEME"></a>
* #### `OPT__REF_MAG_INT_SCHEME` &ensp; (only 2, 3, 4, 6 in ["Supported interpolation schemes"](#INT_TABLE)) &ensp; [4]
    * **Description:**
Interpolation scheme for computing the magnetic field on the newly refined patches.
    * **Restriction:**
For [[--mhd | Installation:-Option-List#--mhd]] only.

<a name="OPT__POT_INT_SCHEME"></a>
* #### `OPT__POT_INT_SCHEME` &ensp; (only 4 & 5 in ["Supported interpolation schemes"](#INT_TABLE)) &ensp; [4]
    * **Description:**
Interpolation scheme for preparing the ghost-zone potential for the Poisson solver.
    * **Restriction:**

<a name="OPT__RHO_INT_SCHEME"></a>
* #### `OPT__RHO_INT_SCHEME` &ensp; (see ["Supported interpolation schemes"](#INT_TABLE)) &ensp; [4]
    * **Description:**
Interpolation scheme for preparing the ghost-zone mass density for the Poisson solver.
    * **Restriction:**

<a name="OPT__GRA_INT_SCHEME"></a>
* #### `OPT__GRA_INT_SCHEME` &ensp; (see ["Supported interpolation schemes"](#INT_TABLE)) &ensp; [4]
    * **Description:**
Interpolation scheme for preparing the ghost-zone potential for computing the
gravitational acceleration.
    * **Restriction:**

<a name="OPT__REF_POT_INT_SCHEME"></a>
* #### `OPT__REF_POT_INT_SCHEME` &ensp; (see ["Supported interpolation schemes"](#INT_TABLE)) &ensp; [4]
    * **Description:**
Interpolation scheme for computing the gravitational potential on the newly refined patches.
    * **Restriction:**

<a name="INT_MONO_COEFF"></a>
* #### `INT_MONO_COEFF` &ensp; (1.0 &#8804; input &#8804; 4.0) &ensp; [2.0]
    * **Description:**
Slope limiter coefficient for ensuring monotonicity in interpolation.
The interpolation results become more diffusive (and presumably also more stable)
when adopting a smaller value.
    * **Restriction:**

<a name="INT_MONO_COEFF_B"></a>
* #### `INT_MONO_COEFF_B` &ensp; (1.0 &#8804; input &#8804; 4.0) &ensp; [2.0]
    * **Description:**
Slope limiter coefficient for ensuring monotonicity when interpolating magnetic field.
The interpolation results become more diffusive (and presumably also more stable)
when adopting a smaller value.
    * **Restriction:**

<a name="MONO_MAX_ITER"></a>
* #### `MONO_MAX_ITER` &ensp; (0=off, &#8805;0=on) &ensp; [10]
    * **Description:**
Maximum number of iterations to reduce [INT_MONO_COEFF](#INT_MONO_COEFF)
when interpolation fails. It improves code stability.
    * **Restriction:**

<a name="INT_OPP_SIGN_0TH_ORDER"></a>
* #### `INT_OPP_SIGN_0TH_ORDER` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Apply 0th-order interpolation if the values to be interpolated change
signs in adjacent cells. This helps avoid introducing unphysically large
velocity in low-density valleys.
    * **Restriction:**


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
