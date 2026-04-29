This page describes the unit system.


## Compilation Options

Related options:
[[--comoving | [Installation]-Option-List#--comoving]] &nbsp;


## Runtime Parameters
[[ [Runtime parameters] Units | [Runtime-Parameters]-Units ]]


## Remarks

### Unit Consistency
When enabling units (i.e., [[OPT__UNIT | [Runtime-Parameters]-Units#OPT__UNIT]])=1), all input values
(e.g., runtime parameters, refinement thresholds, and initial conditions)
must be set consistently with the adopted unit system unless otherwise
specified. All input physical constants (e.g.,
[[NEWTON_G | [Runtime-Parameters]-Gravity#NEWTON_G]],
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
When enabling [[--comoving | [Installation]-Option-List#--comoving]],
the unit system is currently fixed as follows.
- [[UNIT_L | [Runtime-Parameters]-Units#UNIT_L]] = <var>h</var><sup>-1</sup> <var>Mpc</var>,
where <var>h</var>=[[HUBBLE0 | [Runtime-Parameters]-Cosmology#HUBBLE0]] is the dimensionless Hubble constant.
- [[UNIT_T | [Runtime-Parameters]-Units#UNIT_T]] = <var>H</var><sub>0</sub><sup>-1</sup>, where
<var>H</var><sub>0</sub> = <var>100 h km s</var><sup>-1</sup> <var>Mpc</var><sup>-1</sup> is the Hubble constant.
- [[UNIT_D | [Runtime-Parameters]-Units#UNIT_D]] = <var>&rho;</var><sub>m0</sub> =
<var>3&Omega;</var><sub>m0</sub><var>H</var><sub>0</sub><sup>2</sup> / <var>8&pi;G</var>,
where <var>&rho;</var><sub>m0</sub> is the present matter density,
<var>&Omega;</var><sub>m0</sub>=[[OMEGA_M0 | [Runtime-Parameters]-Cosmology#OMEGA_M0]]
is the present matter density parameter, and
<var>G</var> is the gravitational constant.


<br>

## Links
* [[Main page of Runtime Parameters | Runtime-Parameters]]
