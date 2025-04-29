Parameters described on this page:
[FB_MIN_LEVEL](#FB_MIN_LEVEL), &nbsp;
[FB_RSEED](#FB_RSEED), &nbsp;
[FB_SNE](#FB_SNE), &nbsp;
[FB_RESOLVED_SNEII](#FB_RESOLVED_SNEII), &nbsp;
[FB_USER](#FB_USER) &nbsp;
[FB_RESOLVED_SNEII_N_PER_MASS](#FB_RESOLVED_SNEII_N_PER_MASS), &nbsp;
[FB_RESOLVED_SNEII_DELAY_TIME](#FB_RESOLVED_SNEII_DELAY_TIME), &nbsp;
[FB_RESOLVED_SNEII_EJECT_ENGY](#FB_RESOLVED_SNEII_EJECT_ENGY), &nbsp;
[FB_RESOLVED_SNEII_EJECT_MASS](#FB_RESOLVED_SNEII_EJECT_MASS), &nbsp;
[FB_RESOLVED_SNEII_EJECT_METAL](#FB_RESOLVED_SNEII_EJECT_METAL), &nbsp;
[FB_RESOLVED_SNEII_MIN_M_GAS](#FB_RESOLVED_SNEII_MIN_M_GAS), &nbsp;
[FB_RESOLVED_SNEII_RECORD](#FB_RESOLVED_SNEII_RECORD), &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**


<a name="FB_MIN_LEVEL"></a>
* #### `FB_MIN_LEVEL` &ensp; (0 &#8804; input < [[--nlevel | Installation:-Option-List#--nlevel]]; <0 &#8594; set to [[MAX_LEVEL | Runtime Parameters:-Refinement#MAX_LEVEL ]]) &ensp; [0]
    * **Description:**
Minimum AMR level to apply feedback.
    * **Restriction:**

<a name="FB_RSEED"></a>
* #### `FB_RSEED` &ensp; (&#8805;0) &ensp; [456]
    * **Description:**
Random seed used by feedback.
    * **Restriction:**

<a name="FB_SNE"></a>
* #### `FB_SNE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Supernova explosion feedback.
    * **Restriction:**
Not supported yet.

<a name="FB_RESOLVED_SNEII"></a>
* #### `FB_RESOLVED_SNEII` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Resolved Type II supernovae feedback.
When a star particle forms, it would be sampled to contain a SNII progenitor stochastically
with a probability of star particle mass times
number of SNII per stellar masses ([[FB_RESOLVED_SNEII_N_PER_MASS | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII_N_PER_MASS ]]).
The selected particle for supernova will then explode when the age of the star particle reaches
the fixed explosion delay time ([[FB_RESOLVED_SNEII_DELAY_TIME | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII_DELAY_TIME ]]).
When the supernova explodes, it will uniformly inject
the given amount of thermal energy ([[FB_RESOLVED_SNEII_EJECT_ENGY | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII_EJECT_ENGY ]]),
mass ([[FB_RESOLVED_SNEII_EJECT_MASS | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII_EJECT_MASS ]]), and
metal ([[FB_RESOLVED_SNEII_EJECT_METAL | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII_EJECT_METAL ]])
into fluid of the surrounding cells where the particle is located.
The number of cells to apply feedback will increase from one until the enclosed mass exceeds the given
minimum gas mass ([[FB_RESOLVED_SNEII_MIN_M_GAS | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII_MIN_M_GAS ]]), or
the diameter equals `FB_GHOST_SIZE+1`, or the region crosses the coarse-fine boundary.
See sec. 2.6 in [Chia-Yu Hu et al. 2023](https://doi.org/10.3847/1538-4357/accf9e) for reference.
    * **Restriction:**
Must set one extra particle attribute with [[ --par_attribute_flt | Installation:-Option-List#--par_attribute_flt ]].
The star particle mass resolution must be high enough to have at most one supernova explosion per particle.
The grid resolution must be high enough to resolve the Sedov phase of supernova explosion blast wave,
so the kinetic (outward momentum) feedback is not needed.

<a name="FB_USER"></a>
* #### `FB_USER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
User-defined feedback.
See [Add User-defined Feedback](#add-user-defined-feedback) for details.
    * **Restriction:**

<a name="FB_RESOLVED_SNEII_N_PER_MASS"></a>
* #### `FB_RESOLVED_SNEII_N_PER_MASS` &ensp; (&#8805;0.0) &ensp; [1.0e-2]
    * **Description:**
Number of SNeII per stellar mass.
Note that the input value should always be in units of 1/Msun.
    * **Restriction:**
Only for [[FB_RESOLVED_SNEII | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII ]].
Its value times the maximum star particle mass cannot be greater than one.

<a name="FB_RESOLVED_SNEII_DELAY_TIME"></a>
* #### `FB_RESOLVED_SNEII_DELAY_TIME` &ensp; (&#8805;0.0) &ensp; [10.0]
    * **Description:**
Explosion delay time of SNeII after star formation.
Note that the input value should always be in units of Myr.
    * **Restriction:**
Only for [[FB_RESOLVED_SNEII | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII ]].

<a name="FB_RESOLVED_SNEII_EJECT_ENGY"></a>
* #### `FB_RESOLVED_SNEII_EJECT_ENGY` &ensp; (&#8805;0.0) &ensp; [1.0e51]
    * **Description:**
Ejected internal energy of each SNeII explosion.
Note that the input value should always be in units of erg.
    * **Restriction:**
Only for [[FB_RESOLVED_SNEII | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII ]].

<a name="FB_RESOLVED_SNEII_EJECT_MASS"></a>
* #### `FB_RESOLVED_SNEII_EJECT_MASS` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Ejected total mass of each SNeII explosion.
Note that the input value should always be in units of Msun.
The actual amount of the ejected mass will not be greater than the mass of the particle.
    * **Restriction:**
Only for [[FB_RESOLVED_SNEII | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII ]].

<a name="FB_RESOLVED_SNEII_EJECT_METAL"></a>
* #### `FB_RESOLVED_SNEII_EJECT_METAL` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Ejected metal mass of each SNeII explosion.
Note that the input value should always be in units of Msun.
The actual amount of the ejected metal will not be greater than
[[FB_RESOLVED_SNEII_EJECT_MASS | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII_EJECT_MASS ]]
and the metal mass of the particle.
    * **Restriction:**
Only for [[FB_RESOLVED_SNEII | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII ]].

<a name="FB_RESOLVED_SNEII_MIN_M_GAS"></a>
* #### `FB_RESOLVED_SNEII_MIN_M_GAS` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Minimum mass of gas in the environment of SNeII to apply feedback.
If the mass of the cell where the particle is in is sufficient, then the feedback is applied to a single cell.
Otherwise, the diameter of the feedback region will increase from `1` to `FB_GHOST_SIZE+1`
until the enclosed mass is higher than this given threshold.
When the particle is too close to the AMR coarse-fine boundary,
the feedback region will be limited such that no feedback will be applied in the coarse level.
Note that the input value should always be in units of Msun.
    * **Restriction:**
Only for [[FB_RESOLVED_SNEII | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII ]].

<a name="FB_RESOLVED_SNEII_RECORD"></a>
* #### `FB_RESOLVED_SNEII_RECORD` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Record the information of the resolved Type II supernovae feedback.
    * **Restriction:**
Only for [[FB_RESOLVED_SNEII | Runtime Parameters:-Feedback#FB_RESOLVED_SNEII ]].


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Feedback | Feedback]]
