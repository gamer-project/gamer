Parameters described on this page:
[SF_CREATE_STAR_SCHEME](#SF_CREATE_STAR_SCHEME), &nbsp;
[SF_CREATE_STAR_RSEED](#SF_CREATE_STAR_RSEED), &nbsp;
[SF_CREATE_STAR_DET_RANDOM](#SF_CREATE_STAR_DET_RANDOM), &nbsp;
[SF_CREATE_STAR_MIN_LEVEL](#SF_CREATE_STAR_MIN_LEVEL), &nbsp;
[SF_CREATE_STAR_MIN_GAS_DENS](#SF_CREATE_STAR_MIN_GAS_DENS), &nbsp;
[SF_CREATE_STAR_MASS_EFF](#SF_CREATE_STAR_MASS_EFF), &nbsp;
[SF_CREATE_STAR_MIN_STAR_MASS](#SF_CREATE_STAR_MIN_STAR_MASS), &nbsp;
[SF_CREATE_STAR_MAX_STAR_MFRAC](#SF_CREATE_STAR_MAX_STAR_MFRAC), &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="SF_CREATE_STAR_SCHEME"></a>
* #### `SF_CREATE_STAR_SCHEME` &ensp; (0=off, 1=AGORA) &ensp; [0]
    * **Description:**
Star formation schemes.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | Installation:-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_RSEED"></a>
* #### `SF_CREATE_STAR_RSEED` &ensp; (&#8805;0) &ensp; [123]
    * **Description:**
Random seed used by star formation.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | Installation:-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_DET_RANDOM"></a>
* #### `SF_CREATE_STAR_DET_RANDOM` &ensp; (0=off, 1=on; <0=auto) &ensp; [-1]
    * **Description:**
Make random numbers deterministic (i.e., independent of OpenMP and MPI).
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | Installation:-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MIN_LEVEL"></a>
* #### `SF_CREATE_STAR_MIN_LEVEL` &ensp; (&#8805;0; <0 &#8594; MAX_LEVEL) &ensp; [0]
    * **Description:**
Minimum AMR level allowed to form stars.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | Installation:-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MIN_GAS_DENS"></a>
* #### `SF_CREATE_STAR_MIN_GAS_DENS` &ensp; (&#8805;0.0) &ensp; [1.0e1]
    * **Description:**
Minimum gas density allowed to form stars.
Note that the input value should always be in units of HI count/cm^3.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | Installation:-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MASS_EFF"></a>
* #### `SF_CREATE_STAR_MASS_EFF` &ensp; (>0.0) &ensp; [1.0e-2]
    * **Description:**
Gas-to-star mass conversion efficiency.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | Installation:-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MIN_STAR_MASS"></a>
* #### `SF_CREATE_STAR_MIN_STAR_MASS` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Minimum star particle mass for the stochastical star formation.
Note that the input value should always be in units of Msun.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | Installation:-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MAX_STAR_MFRAC"></a>
* #### `SF_CREATE_STAR_MAX_STAR_MFRAC` &ensp; (>0.0) &ensp; [0.5]
    * **Description:**
Maximum gas mass fraction allowed to convert to stars per substep.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | Installation:-Option-List#--star_formation]].



## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Star Formation | Star-Formation]]
