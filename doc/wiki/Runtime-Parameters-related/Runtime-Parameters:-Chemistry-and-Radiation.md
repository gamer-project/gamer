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
[GRACKLE_CLOUDY_TABLE](#GRACKLE_CLOUDY_TABLE) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="GRACKLE_ACTIVATE"></a>
* #### `GRACKLE_ACTIVATE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Enable the [GRACKLE chemistry and cooling library](http://grackle.readthedocs.io/en/latest/index.html).
    * **Restriction:**
Only applicable when enabling the compilation option
[[--grackle | Installation:-Option-List#--grackle]].

<a name="GRACKLE_VERBOSE"></a>
* #### `GRACKLE_VERBOSE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Map to the "grackle_verbose" option in GRACKLE.
    * **Restriction:**

<a name="GRACKLE_COOLING"></a>
* #### `GRACKLE_COOLING` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Map to the "with_radiative_cooling" runtime parameter in GRACKLE.
    * **Restriction:**

<a name="GRACKLE_PRIMORDIAL"></a>
* #### `GRACKLE_PRIMORDIAL` &ensp; (0=Cloudy, 1=6-species, 2=9-species, 3=12-species) &ensp; [0]
    * **Description:**
Map to the "primordial_chemistry" runtime parameter in GRACKLE.
One must increase
[[--passive | Installation:-Option-List#--passive]]
by 3, 6, or 9 for GRACKLE_PRIMORDIAL=1, 2, or 3, respectively.
    * **Restriction:**

<a name="GRACKLE_METAL"></a>
* #### `GRACKLE_METAL` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Map to the "metal_cooling" runtime parameter in GRACKLE. One must increase
[[--passive | Installation:-Option-List#--passive]]
by 1 and initialize the field `Metal` using the field index `Idx_Metal` properly.
    * **Restriction:**

<a name="GRACKLE_UV"></a>
* #### `GRACKLE_UV` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Map to the "UVbackground" runtime parameter in GRACKLE.
    * **Restriction:**

<a name="GRACKLE_CMB_FLOOR"></a>
* #### `GRACKLE_CMB_FLOOR` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Map to the "cmb_temperature_floor" runtime parameter in GRACKLE.
    * **Restriction:**

<a name="GRACKLE_PE_HEATING"></a>
* #### `GRACKLE_PE_HEATING` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Map to the "photoelectric_heating" runtime parameter in GRACKLE.
    * **Restriction:**

<a name="GRACKLE_PE_HEATING_RATE"></a>
* #### `GRACKLE_PE_HEATING_RATE` &ensp; (&#8805;0.0) &ensp; [8.5e-26]
    * **Description:**
Map to the "photoelectric_heating_rate" runtime parameter in GRACKLE.
Note that the input value should always be in units of
<var>erg</var>&#8287;<var>cm</var><sup>-3</sup>&#8287;<var>s</var><sup>-1</sup>.
    * **Restriction:**

<a name="GRACKLE_CLOUDY_TABLE"></a>
* #### `GRACKLE_CLOUDY_TABLE` &ensp; (string) &ensp; [none]
    * **Description:**
Map to the "grackle_data_file" runtime parameter in GRACKLE.
    * **Restriction:**


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Chemistry and Radiation | Chemistry-and-Radiation]]
