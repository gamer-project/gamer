Parameters described on this page:
[OPT__UNIT](#OPT__UNIT), &nbsp;
[UNIT_L](#UNIT_L), &nbsp;
[UNIT_M](#UNIT_M), &nbsp;
[UNIT_T](#UNIT_T), &nbsp;
[UNIT_V](#UNIT_V), &nbsp;
[UNIT_D](#UNIT_D) &nbsp;

Other related parameters:
[[NEWTON_G| [Runtime-Parameters]-Gravity#NEWTON_G]], &nbsp;
[[ELBDM_MASS | Wave-Dark-Matter#ELBDM_MASS]], &nbsp;
[[ELBDM_PLANCK_CONST | Wave-Dark-Matter#ELBDM_PLANCK_CONST]] &nbsp;

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="OPT__UNIT"></a>
* #### `OPT__UNIT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Enable code units. One must also set exactly *three* units among
[UNIT_L](#UNIT_L), [UNIT_M](#UNIT_M), [UNIT_T](#UNIT_T), [UNIT_V](#UNIT_V), [UNIT_D](#UNIT_D)
(unless the compilation option [[--comoving | [Installation]-Option-List#--comoving]]
is adopted; see [[Units in Cosmological Simulations | Units#units-in-cosmological-simulations]] for details).
See also [[Unit Consistency | Units#unit-consistency]].
    * **Restriction:**
It will be turned on automatically when enabling the compilation option
[[--comoving | [Installation]-Option-List#--comoving]].

<a name="UNIT_L"></a>
* #### `UNIT_L` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Length unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | [Installation]-Option-List#--comoving]].
See [[Units in Cosmological Simulations | Units#units-in-cosmological-simulations]]
for details.

<a name="UNIT_M"></a>
* #### `UNIT_M` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Mass unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | [Installation]-Option-List#--comoving]].
See [[Units in Cosmological Simulations | Units#units-in-cosmological-simulations]]
for details.

<a name="UNIT_T"></a>
* #### `UNIT_T` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Time unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | [Installation]-Option-List#--comoving]].
See [[Units in Cosmological Simulations | Units#units-in-cosmological-simulations]]
for details.

<a name="UNIT_V"></a>
* #### `UNIT_V` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Velocity unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | [Installation]-Option-List#--comoving]].
See [[Units in Cosmological Simulations | Units#units-in-cosmological-simulations]]
for details.

<a name="UNIT_D"></a>
* #### `UNIT_D` &ensp; (>0.0; &#8804;0.0 &#8594; set by other units) &ensp; [none]
    * **Description:**
Mass density unit in CGS.
    * **Restriction:**
It will be overwritten when enabling the compilation option
[[--comoving | [Installation]-Option-List#--comoving]].
See [[Units in Cosmological Simulations | Units#units-in-cosmological-simulations]]
for details.


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime-Parameters]]
* [[Main page of Units | Units]]
