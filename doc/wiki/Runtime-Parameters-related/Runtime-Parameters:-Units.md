
## Compilation Options

Related options:
[[--comoving | Installation:-Option-List#--comoving]] &nbsp;

## Runtime Parameters

Parameters described on this page:
[OPT__UNIT](#OPT__UNIT), &nbsp;
[UNIT_L](#UNIT_L), &nbsp;
[UNIT_M](#UNIT_M), &nbsp;
[UNIT_T](#UNIT_T), &nbsp;
[UNIT_V](#UNIT_V), &nbsp;
[UNIT_D](#UNIT_D) &nbsp;

Other related parameters:
[[NEWTON_G| Runtime-Parameters:-Gravity#NEWTON_G]], &nbsp;
[[ELBDM_MASS | Wave-Dark-Matter#ELBDM_MASS]], &nbsp;
[[ELBDM_PLANCK_CONST | Wave-Dark-Matter#ELBDM_PLANCK_CONST]] &nbsp;

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="OPT__UNIT"></a>
* #### `OPT__UNIT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Enable code units. One must also set exactly *three* units among
[UNIT_L](#UNIT_L), [UNIT_M](#UNIT_M), [UNIT_T](#UNIT_T), [UNIT_V](#UNIT_V), [UNIT_D](#UNIT_D)
(unless the compilation option [[--comoving | Installation:-Option-List#--comoving]]
is adopted; see [Units in Cosmological Simulations](#units-in-cosmological-simulations) for details).
See also [Unit Consistency](#unit-consistency).
    * **Restriction:**
It will be turned on automatically when enabling the compilation option
[[--comoving | Installation:-Option-List#--comoving]].

<a name="UNIT_L"></a>
* #### `UNIT_L` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Length unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | Installation:-Option-List#--comoving]].
See [Units in Cosmological Simulations](#units-in-cosmological-simulations)
for details.

<a name="UNIT_M"></a>
* #### `UNIT_M` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Mass unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | Installation:-Option-List#--comoving]].
See [Units in Cosmological Simulations](#units-in-cosmological-simulations)
for details.

<a name="UNIT_T"></a>
* #### `UNIT_T` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Time unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | Installation:-Option-List#--comoving]].
See [Units in Cosmological Simulations](#units-in-cosmological-simulations)
for details.

<a name="UNIT_V"></a>
* #### `UNIT_V` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Velocity unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | Installation:-Option-List#--comoving]].
See [Units in Cosmological Simulations](#units-in-cosmological-simulations)
for details.

<a name="UNIT_D"></a>
* #### `UNIT_D` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Mass density unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | Installation:-Option-List#--comoving]].
See [Units in Cosmological Simulations](#units-in-cosmological-simulations)
for details.


## Remarks

### Unit Consistency
When enabling units (i.e., [OPT__UNIT](#OPT__UNIT)=1), all input values
(e.g., runtime parameters, refinement thresholds, and initial conditions)
must be set consistently with the adopted unit system unless otherwise
specified. All input physical constants (e.g.,
[[NEWTON_G| Runtime-Parameters:-Gravity#NEWTON_G]],
[[ELBDM_MASS | Wave-Dark-Matter#ELBDM_MASS]], and
[[ELBDM_PLANCK_CONST | Wave-Dark-Matter#ELBDM_PLANCK_CONST]])
will be reset automatically to conform to the adopted unit system.

The code does not distinguish external and internal units;
In other words, all internal variables in the code should also conform to
the same unit system (except for the predefined
[physical constants](#physical-constants)).

### Physical Constants
Some physical constants in CGS units are predefined in the header
`include/PhysicalConstant.h` with the prefix `Const_`. They can be
useful for unit conversion.
For example, to convert particle velocity `vel_code` from code units
to km/s, one can use `vel_kms = vel_code*UNIT_V*Const_s/Const_km;`.

### Units in Cosmological Simulations
When enabling [[--comoving | Installation:-Option-List#--comoving]],
the unit system is currently fixed as follows.
- [UNIT_L](#UNIT_L) = <var>h</var><sup>-1</sup> <var>Mpc</var>,
where <var>h</var>=[[HUBBLE0 | Runtime-Parameters:-Cosmology#HUBBLE0]] is the dimensionless Hubble constant.
- [UNIT_T](#UNIT_T) = <var>H</var><sub>0</sub><sup>-1</sup>, where
<var>H</var><sub>0</sub> = <var>100 h km s</var><sup>-1</sup> <var>Mpc</var><sup>-1</sup> is the Hubble constant.
- [UNIT_D](#UNIT_D) = <var>&rho;</var><sub>m0</sub> =
<var>3&Omega;</var><sub>m0</sub><var>H</var><sub>0</sub><sup>2</sup> / <var>8&pi;G</var>,
where <var>&rho;</var><sub>m0</sub> is the present matter density,
<var>&Omega;</var><sub>m0</sub>=[[OMEGA_M0 | Runtime-Parameters:-Cosmology#OMEGA_M0]]
is the present matter density parameter, and
<var>G</var> is the gravitational constant.


<br>

## Links
* [[Main page of Runtime Parameters | Runtime-Parameters]]
