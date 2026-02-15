Parameters described on this page:
[GRACKLE_ACTIVATE](#GRACKLE_ACTIVATE), &nbsp;
[GRACKLE_VERBOSE](#GRACKLE_VERBOSE), &nbsp;
[GRACKLE_COOLING](#GRACKLE_COOLING), &nbsp;
[GRACKLE_PRIMORDIAL](#GRACKLE_PRIMORDIAL), &nbsp;
[GRACKLE_METAL](#GRACKLE_METAL), &nbsp;
[GRACKLE_UV](#GRACKLE_UV), &nbsp;
[GRACKLE_CMB_FLOOR](#GRACKLE_CMB_FLOOR), &nbsp;
[GRACKLE_PE_HEATING](#GRACKLE_PE_HEATING), &nbsp;
[GRACKLE_PE_HEATING_RATE](#GRACKLE_PE_HEATING_RATE), &nbsp;
[GRACKLE_CLOUDY_TABLE](#GRACKLE_CLOUDY_TABLE), &nbsp;
[GRACKLE_THREE_BODY_RATE](#GRACKLE_THREE_BODY_RATE), &nbsp;
[GRACKLE_CIE_COOLING](#GRACKLE_CIE_COOLING), &nbsp;
[GRACKLE_H2_OPA_APPROX](#GRACKLE_H2_OPA_APPROX), &nbsp;
[GRACKLE_USE_V_HEATING_RATE](#GRACKLE_USE_V_HEATING_RATE), &nbsp;
[GRACKLE_USE_S_HEATING_RATE](#GRACKLE_USE_S_HEATING_RATE), &nbsp;
[GRACKLE_USE_TEMP_FLOOR](#GRACKLE_USE_TEMP_FLOOR), &nbsp;
[GRACKLE_TEMP_FLOOR_SCALAR](#GRACKLE_TEMP_FLOOR_SCALAR), &nbsp;
[GRACKLE_HYDROGEN_MFRAC](#GRACKLE_HYDROGEN_MFRAC), &nbsp;
[OPT__UNFREEZE_GRACKLE](#OPT__UNFREEZE_GRACKLE) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="GRACKLE_ACTIVATE"></a>
* #### `GRACKLE_ACTIVATE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Enable the [GRACKLE chemistry and cooling library](http://grackle.readthedocs.io/en/latest/index.html).
    * **Restriction:**
Only applicable when enabling the compilation option
[[--grackle | [Installation]-Option-List#--grackle]].

<a name="GRACKLE_VERBOSE"></a>
* #### `GRACKLE_VERBOSE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Map to the ["grackle_verbose" option in GRACKLE](https://grackle.readthedocs.io/en/latest/Interaction.html#enabling-output).
    * **Restriction:**

<a name="GRACKLE_COOLING"></a>
* #### `GRACKLE_COOLING` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Map to the ["with_radiative_cooling" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.with_radiative_cooling).
    * **Restriction:**

<a name="GRACKLE_PRIMORDIAL"></a>
* #### `GRACKLE_PRIMORDIAL` &ensp; (0=Cloudy, 1=6-species, 2=9-species, 3=12-species) &ensp; [0]
    * **Description:**
Map to the ["primordial_chemistry" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.primordial_chemistry).
One must increase
[[--passive | [Installation]-Option-List#--passive]]
by 6, 9, or 12 for GRACKLE_PRIMORDIAL=1, 2, or 3, respectively.
    * **Restriction:**

<a name="GRACKLE_METAL"></a>
* #### `GRACKLE_METAL` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Map to the ["metal_cooling" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.metal_cooling). One must increase
[[--passive | [Installation]-Option-List#--passive]]
by 1 and initialize the field `Metal` using the field index `Idx_Metal` properly.
    * **Restriction:**

<a name="GRACKLE_UV"></a>
* #### `GRACKLE_UV` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Map to the ["UVbackground" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.UVbackground).
    * **Restriction:**

<a name="GRACKLE_CMB_FLOOR"></a>
* #### `GRACKLE_CMB_FLOOR` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Map to the ["cmb_temperature_floor" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.cmb_temperature_floor).
    * **Restriction:**

<a name="GRACKLE_PE_HEATING"></a>
* #### `GRACKLE_PE_HEATING` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Map to the ["photoelectric_heating" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.photoelectric_heating).
    * **Restriction:**

<a name="GRACKLE_PE_HEATING_RATE"></a>
* #### `GRACKLE_PE_HEATING_RATE` &ensp; (&#8805;0.0) &ensp; [8.5e-26]
    * **Description:**
Map to the ["photoelectric_heating_rate" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.photoelectric_heating_rate).
Note that the input value should always be in units of
<var>erg</var>&#8287;<var>cm</var><sup>-3</sup>&#8287;<var>s</var><sup>-1</sup>&#8287;<var>n</var><sub>H</sub><sup>-1</sup>.
    * **Restriction:**

<a name="GRACKLE_CLOUDY_TABLE"></a>
* #### `GRACKLE_CLOUDY_TABLE` &ensp; (string) &ensp; [none]
    * **Description:**
Map to the ["grackle_data_file" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.grackle_data_file).
    * **Restriction:**

<a name="GRACKLE_THREE_BODY_RATE"></a>
* #### `GRACKLE_THREE_BODY_RATE` &ensp; (0=Abel+02, 1=Palla+83, 2=Cohen+83, 3=Flower+07, 4=Glover+08, 5=Forrey+13) &ensp; [0]
    * **Description:**
Map to the ["three_body_rate" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.three_body_rate).
    * **Restriction:**

<a name="GRACKLE_CIE_COOLING"></a>
* #### `GRACKLE_CIE_COOLING` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Map to the ["cie_cooling" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.cie_cooling).
    * **Restriction:**

<a name="GRACKLE_H2_OPA_APPROX"></a>
* #### `GRACKLE_H2_OPA_APPROX` &ensp; (0=off, 1=Ripomonti+04) &ensp; [0]
    * **Description:**
Map to the ["h2_optical_depth_approximation" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.h2_optical_depth_approximation).
    * **Restriction:**

<a name="GRACKLE_USE_V_HEATING_RATE"></a>
* #### `GRACKLE_USE_V_HEATING_RATE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Map to the ["use_volumetric_heating_rate" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.use_volumetric_heating_rate).
    * **Restriction:**

<a name="GRACKLE_USE_S_HEATING_RATE"></a>
* #### `GRACKLE_USE_S_HEATING_RATE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Map to the ["use_specific_heating_rate" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.use_specific_heating_rate).
    * **Restriction:**

<a name="GRACKLE_USE_TEMP_FLOOR"></a>
* #### `GRACKLE_USE_TEMP_FLOOR` &ensp; (0=off, 1=single value, 2=array of values) &ensp; [0]
    * **Description:**
Map to the ["use_temperature_floor" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.use_temperature_floor).
For `GRACKLE_USE_TEMP_FLOOR == 1`, set [GRACKLE_TEMP_FLOOR_SCALAR](#GRACKLE_TEMP_FLOOR_SCALAR) as the single value.
For `GRACKLE_USE_TEMP_FLOOR == 2`, define your function in the test problem and link it to `Grackle_tempFloor_User_Ptr` to set the array of values.
    * **Restriction:**

<a name="GRACKLE_TEMP_FLOOR_SCALAR"></a>
* #### `GRACKLE_TEMP_FLOOR_SCALAR` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Map to the ["temperature_floor_scalar" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.temperature_floor_scalar).
    * **Restriction:**
It is used only when [GRACKLE_USE_TEMP_FLOOR](#GRACKLE_USE_TEMP_FLOOR)==1.

<a name="GRACKLE_HYDROGEN_MFRAC"></a>
* #### `GRACKLE_PE_HEATING_RATE` &ensp; (&#8805;0.0) &ensp; [0.76]
    * **Description:**
Map to the ["HydrogenFractionByMass" runtime parameter in GRACKLE](https://grackle.readthedocs.io/en/latest/Parameters.html#c.HydrogenFractionByMass).
Note that the input value will only set to Grackle when it is in non-equilibrium mode [GRACKLE_PRIMORDIAL](#GRACKLE_PRIMORDIAL)>0,
because the tables for tabulated mode were created assuming hydrogen mass fraction of about 0.716.
See the Grackle document for the details.
When in tabulated mode ([GRACKLE_PRIMORDIAL](#GRACKLE_PRIMORDIAL)==0), this value may still be used elsewhere in GAMER and
it may be different from HydrogenFractionByMass inside Grackle.
    * **Restriction:**

<a name="OPT__UNFREEZE_GRACKLE"></a>
* #### `OPT__UNFREEZE_GRACKLE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Allow the evolution by Grackle solver when the fluid is frozen
(for [[OPT__FREEZE_FLUID | [Runtime-Parameters]-Hydro#OPT__FREEZE_FLUID]]==1 only).
    * **Restriction:**


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Chemistry and Radiation | Chemistry-and-Radiation]]
