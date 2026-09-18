Parameters described on this page:
[SRC_DELEPTONIZATION](#SRC_DELEPTONIZATION), &nbsp;
[SRC_TURBULENCE](#SRC_TURBULENCE), &nbsp;
[SRC_TURB_VEL](#SRC_TURB_VEL), &nbsp;
[SRC_TURB_AMPL_FACTOR](#SRC_TURB_AMPL_FACTOR), &nbsp;
[SRC_TURB_KDRIV](#SRC_TURB_KDRIV), &nbsp;
[SRC_TURB_KMIN](#SRC_TURB_KMIN), &nbsp;
[SRC_TURB_KMAX](#SRC_TURB_KMAX), &nbsp;
[SRC_TURB_ZETA](#SRC_TURB_ZETA), &nbsp;
[SRC_TURB_SPEC_FORM](#SRC_TURB_SPEC_FORM), &nbsp;
[SRC_TURB_POW](#SRC_TURB_POW), &nbsp;
[SRC_TURB_RSEED_INIT](#SRC_TURB_RSEED_INIT), &nbsp;
[SRC_TURB_UPDATE_STEP](#SRC_TURB_UPDATE_STEP), &nbsp;
[SRC_TURB_TABLE_SIZE](#SRC_TURB_TABLE_SIZE), &nbsp;
[SRC_TURB_RESET](#SRC_TURB_RESET), &nbsp;
[SRC_USER](#SRC_USER) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="SRC_DELEPTONIZATION"></a>
* #### `SRC_DELEPTONIZATION` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Deleptonization (for simulations of stellar core collapse).
    * **Restriction:**
Only applicable when enabling the compilation option
[[--model | [Installation]-Option-List#--model]]=HYDRO.

<a name="SRC_TURBULENCE"></a>
* #### `SRC_TURBULENCE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Turbulence acceleration source term.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--model | [Installation]-Option-List#--model]]=HYDRO.
Does not work with [[--comoving | [Installation]-Option-List#--comoving]] enabled.
Currently only support cubic simulation domain.

<a name="SRC_TURB_VEL"></a>
* #### `SRC_TURB_VEL` &ensp; (>0.0) &ensp; [0.2]
    * **Description:**
Target turbulence velocity dispersion in code unit.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_AMPL_FACTOR"></a>
* #### `SRC_TURB_AMPL_FACTOR` &ensp; (>0.0) &ensp; [1.0]
    * **Description:**
Amplitude multiplier for turbulence acceleration field.
Adjust this factor if velocity dispersion measured is different from target velocity dispersion.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_KDRIV"></a>
* #### `SRC_TURB_KDRIV` &ensp; (>0.0) &ensp; [2.0]
    * **Description:**
Turbulence driving wave number in units of 2π/BOX_SIZE, determine the correlation time τ=BOX_SIZE/(SRC_TURB_KDRIV*SRC_TURB_VEL)
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_KMIN"></a>
* #### `SRC_TURB_KMIN` &ensp; (&#8805;1.0) &ensp; [1.0]
    * **Description:**
Turbulence minimum wave number in units of 2π/BOX_SIZE.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_KMAX"></a>
* #### `SRC_TURB_KMAX` &ensp; (&#8805;1.0) &ensp; [3.0]
    * **Description:**
Turbulence maximum wave number in units of 2π/BOX_SIZE.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_ZETA"></a>
* #### `SRC_TURB_ZETA` &ensp; (1.0&#8805;input&#8805;0.0) &ensp; [1.0]
    * **Description:**
Turbulence solenoidal weighting (0.0=fully compressive, 1.0=fully solenoidal).
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_SPEC_FORM"></a>
* #### `SRC_TURB_SPEC_FORM` &ensp; (0=constant, 1=parabolic, 2=power law) &ensp; [1]
    * **Description:**
Spectral form of the turbulence driving amplitude.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_POW"></a>
* #### `SRC_TURB_POW` &ensp; (real value) &ensp; [-5/3]
    * **Description:**
Power law for energy power spectrum.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1, [SRC_TURB_SPEC_FORM](#SRC_TURB_SPEC_FORM)=2.

<a name="SRC_TURB_RSEED_INIT"></a>
* #### `SRC_TURB_RSEED_INIT` &ensp; (&#8805;0) &ensp; [123]
    * **Description:**
Initial random seed for turbulence, useless when restart without reset.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_UPDATE_STEP"></a>
* #### `SRC_TURB_UPDATE_STEP` &ensp; (&#8805;1) &ensp; [10]
    * **Description:**
Update turbulence pattern every dt=(correlation time/SRC_TURB_UPDATE_STEP).
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_TABLE_SIZE"></a>
* #### `SRC_TURB_TABLE_SIZE` &ensp; (&#8805;1, must be a power of 2) &ensp; [128]
    * **Description:**
Resolution of turbulence table, must be a power of 2.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_TURB_RESET"></a>
* #### `SRC_TURB_RESET` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Reset turbulence field when restart.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="SRC_USER"></a>
* #### `SRC_USER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
User-defined source terms.
See [[Add User-defined Source Terms | Source-Terms#add-user-defined-source-terms]] for details.
    * **Restriction:**


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Source Terms | Source-Terms]]
