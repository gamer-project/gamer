
Parameters described on this page:
[BOX_SIZE](#BOX_SIZE), &nbsp;
[NX0_TOT_X](#NX0_TOT_X), &nbsp;
[NX0_TOT_Y](#NX0_TOT_Y), &nbsp;
[NX0_TOT_Z](#NX0_TOT_Z), &nbsp;
[END_T](#END_T), &nbsp;
[END_STEP](#END_STEP), &nbsp;
[TESTPROB_ID](#TESTPROB_ID) &nbsp;

Other related parameters: none

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="BOX_SIZE"></a>
* #### `BOX_SIZE` &ensp; (>0.0) &ensp; [none]
    * **Description:**
Simulation box length along the longest dimension of
([NX0_TOT_X](#NX0_TOT_X), [NX0_TOT_Y](#NX0_TOT_Y), [NX0_TOT_Z](#NX0_TOT_Z)).
It must conform to the adopted unit system
(see [[Unit Consistency | Runtime Parameters:-Units#unit-consistency]]).
    * **Restriction:**

<a name="NX0_TOT_X"></a>
* #### `NX0_TOT_X` &ensp; (&#8805;16) &ensp; [none]
    * **Description:**
Number of root-level cells along x.
    * **Restriction:**
Must be a multiple of 16 (i.e., two patches).

<a name="NX0_TOT_Y"></a>
* #### `NX0_TOT_Y` &ensp; (&#8805;16) &ensp; [none]
    * **Description:**
Number of root-level cells along y.
    * **Restriction:**
Must be a multiple of 16 (i.e., two patches).

<a name="NX0_TOT_Z"></a>
* #### `NX0_TOT_Z` &ensp; (&#8805;16) &ensp; [none]
    * **Description:**
Number of root-level cells along z.

    * **Restriction:**
Must be a multiple of 16 (i.e., two patches).

<a name="END_T"></a>
* #### `END_T` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [depend]
    * **Description:**
Simulation end time. `END_T<0.0` is allowed only during restart
(i.e., [[OPT__INIT | Runtime-Parameters:-Initial-Conditions#OPT__INIT]]=2, for which it
will be reset to the value stored in the restart file) or if it will be
reset by the adopted test problem. It must conform to the adopted unit
system (see [[Unit Consistency | Runtime Parameters:-Units#unit-consistency]]).
    * **Restriction:**

<a name="END_STEP"></a>
* #### `END_STEP` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [depend]
    * **Description:**
Simulation end root-level step. `END_STEP<0` is allowed only during restart
(i.e., [[OPT__INIT | Runtime-Parameters:-Initial-Conditions#OPT__INIT]]=2, for which it
will be reset to the value stored in the restart file) or if it will be
reset by the adopted test problem. For `END_STEP=0`, the program will
still construct and output the initial condition before termination.
    * **Restriction:**

<a name="TESTPROB_ID"></a>
* #### `TESTPROB_ID` &ensp; (&#8805;0) &ensp; [0]
    * **Description:**
Test problem ID. See [[Test Problems]] for details. It currently supports
        * 0: none
        * 1: blast wave [HYDRO]
        * 2: acoustic wave [HYDRO]
        * 4: cluster merger [HYDRO + GRAVITY + PARTICLE]
        * 5: AGORA isolated disk galaxy [HYDRO + GRAVITY + PARTICLE + SUPPORT_GRACKLE + STAR_FORMATION]
        * 6: caustic wave [HYDRO]
        * 7: spherical collapse [HYDRO + GRAVITY + COMOVING]
        * 8: Kelvin Helmholtz instability [HYDRO]
        * 9: Riemann problems [HYDRO]
        * ...
    * **Restriction:**


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
