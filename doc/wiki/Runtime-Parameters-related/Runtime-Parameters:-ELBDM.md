Parameters described on this page:
[ELBDM_MASS](#ELBDM_MASS), &nbsp;
[ELBDM_PLANCK_CONST](#ELBDM_PLANCK_CONST), &nbsp;
[ELBDM_LAMBDA](#ELBDM_LAMBDA), &nbsp;
[ELBDM_TAYLOR3_COEFF](#ELBDM_TAYLOR3_COEFF), &nbsp;
[ELBDM_TAYLOR3_AUTO](#ELBDM_TAYLOR3_AUTO), &nbsp;
[ELBDM_REMOVE_MOTION_CM](#ELBDM_REMOVE_MOTION_CM), &nbsp;
[ELBDM_BASE_SPECTRAL](#ELBDM_BASE_SPECTRAL), &nbsp;
[ELBDM_MATCH_PHASE](#ELBDM_MATCH_PHASE), &nbsp;
[ELBDM_FIRST_WAVE_LEVEL](#ELBDM_FIRST_WAVE_LEVEL), &nbsp;
[OPT__RES_PHASE](#OPT__RES_PHASE), &nbsp;
[SPEC_INT_TABLE_PATH](#SPEC_INT_TABLE_PATH), &nbsp;
[SPEC_INT_XY_INSTEAD_DEPHA](#SPEC_INT_XY_INSTEAD_DEPHA), &nbsp;
[SPEC_INT_VORTEX_THRESHOLD](#SPEC_INT_VORTEX_THRESHOLD), &nbsp;
[SPEC_INT_GHOST_BOUNDARY](#SPEC_INT_GHOST_BOUNDARY) &nbsp;

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="ELBDM_MASS"></a>
* #### `ELBDM_MASS` &ensp; (>0) &ensp; [none]
    * **Description:** 
Particle mass in ev/c^2.
(Input unit is fixed even when OPT__UNIT or COMOVING is on)
    * **Restriction:**

<a name="ELBDM_PLANCK_CONST"></a>
* #### `ELBDM_PLANCK_CONST` &ensp; (>0) &ensp; [conform to the unit system set by [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]] or [[--comoving | Installation:-Option-List#--comoving]]]
    * **Description:** 
Reduced planck constant in g.cm^2/s^2.
    * **Restriction:**
It will be overwritten by the default value when [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]] or [[--comoving | Installation:-Option-List#--comoving]]
is on; no default when [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]] and
[[--comoving | Installation:-Option-List#--comoving]] is off.

<a name="ELBDM_LAMBDA"></a>
* #### `ELBDM_LAMBDA` &ensp; (none) &ensp; [1]
    * **Description:**
Quartic self-interaction coefficient in ELBDM.
    * **Restriction:**

<a name="ELBDM_TAYLOR3_COEFF"></a>
* #### `ELBDM_TAYLOR3_COEFF` &ensp; (&#8805;0.125) &ensp; [1.0/6.0]
    * **Description:**
Coefficient for the 3rd-order Taylor expansion of the wave function.
Values below 0.125 are always unstable.
Values &#8804; 1/6 become unstable if
[[DT__FLUID | Runtime-Parameters:-Timestep#DT__FLUID]] >
$\sqrt{3}\pi/8$ or $\sqrt{27}\pi/32$ (when [[--laplacian_four | Installation:-Option-List#--laplacian_four]] is enabled).
    * **Restriction:**
Only applicable when the compilation option
[[ --wave_scheme | Installation:-Option-List#--wave_scheme ]] = `FD`.
Ignored if [ELBDM_TAYLOR3_AUTO](#ELBDM_TAYLOR3_AUTO) is enable.

* #### `ELBDM_TAYLOR3_AUTO` &ensp; (none) &ensp; [0]
    * **Description:**
If this parameter is set to 1, the code will automatically determine the coefficient
[ELBDM_TAYLOR3_COEFF](#ELBDM_TAYLOR3_COEFF) to minimize the amplitude error
for the smallest wavelength.
    * **Restriction:**
Useless if [[ OPT__FREEZE_FLUID | Runtime-Parameters:-Hydro#OPT__FREEZE_FLUID]] is on.

<a name="ELBDM_REMOVE_MOTION_CM"></a>
* #### `ELBDM_REMOVE_MOTION_CM` &ensp; (0=none, 1=init, 2=every step) &ensp; [0]
    * **Description:**
Remove the motion of center-of-mass.
    * **Restriction:**
Only applicable when enabled
[[ OPT__CK_CONSERVATION | Runtime-Parameters:-Miscellaneous#OPT__CK_CONSERVATION ]].
Not supported when
[[ --bitwise_reproducibility | Installation:-Option-List#--bitwise_reproducibility ]] = true.

<a name="ELBDM_BASE_SPECTRAL"></a>
* #### `ELBDM_BASE_SPECTRAL` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Adopt the spectral method to evolve base-level wave function.
    * **Restriction:**
Requires [[--fftw | Installation:-Option-List#--fftw]] = FFTW2/FFTW
and periodic boundary conditions for all directions:
[[OPT__BC_FLU | Runtime-Parameters:-Hydro#OPT__BC_FLU_XM]] = 1.

<a name="ELBDM_MATCH_PHASE"></a>
* #### `ELBDM_MATCH_PHASE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Match child phases with father phases during data restriction.
    * **Restriction:**
Only applicable when enabling the compilation option
[[ --elbdm_scheme | Installation:-Option-List#--elbdm_scheme]] = `HYBRID`.
Requires [[ OPT__UM_IC_LEVEL | Runtime-Parameters:-Initial-Conditions#OPT__UM_IC_LEVEL ]]
&#8805; [ELBDM_FIRST_WAVE_LEVEL](#ELBDM_FIRST_WAVE_LEVEL).

<a name="ELBDM_FIRST_WAVE_LEVEL"></a>
* #### `ELBDM_FIRST_WAVE_LEVEL` &ensp; (1 &#8804; input &#8804; [[ MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]) &ensp; [none]
    * **Description:**
Level at which to switch to the wave solver.
    * **Restriction:**
Only applicable when enabling the compilation option
[[ --elbdm_scheme | Installation:-Option-List#--elbdm_scheme]] = `HYBRID`.

<a name="OPT__RES_PHASE"></a>
* #### `OPT__RES_PHASE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Restriction on phase. (i.e., the averages of fine-grid data equal the coarse-grid data)
    * **Restriction:**

<a name="SPEC_INT_TABLE_PATH"></a>
* #### `SPEC_INT_TABLE_PATH` &ensp; (none) &ensp; [none]
    * **Description:**
Path to the table of the spectral interpolation.
See [[ ELBDM Spectral Interpolation | ELBDM-Spectral-Interpolation]] for details.
Table download script is available at
`example/test_problem/ELBDM/LSS_Hybrid/download_spectral_interpolation_tables.sh`.
    * **Restriction:**
Only applicable when the enabling compilation option
[[ --spectral_interpolation | Installation:-Option-List#--spectral_interpolation]]
and [[ Interpolation schemes | Runtime-Parameters:-Interpolation##INT_TABLE]] = 8.

<a name="SPEC_INT_XY_INSTEAD_DEPHA"></a>
* #### `SPEC_INT_XY_INSTEAD_DEPHA` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Interpolate x and y (real and imaginary parts in current implementation) around vortices
instead of density and phase for the spectral interpolation,
which has the advantage of being well-defined across vortices
    * **Restriction:**
Only applicable when the enabling compilation option
[[ --spectral_interpolation | Installation:-Option-List#--spectral_interpolation]]
and [[ Interpolation schemes | Runtime-Parameters:-Interpolation##INT_TABLE]] = 8.

<a name="SPEC_INT_VORTEX_THRESHOLD"></a>
* #### `SPEC_INT_VORTEX_THRESHOLD` &ensp; (&#8805;0) &ensp; [0.1]
    * **Description:**
Vortex detection threshold for [SPEC_INT_XY_INSTEAD_DEPHA](#SPEC_INT_XY_INSTEAD_DEPHA),
triggered when Lap(S) * dx**2 > threshold, indicating a significant phase jump.
    * **Restriction:**
Only applicable when the enabling compilation option
[[ --spectral_interpolation | Installation:-Option-List#--spectral_interpolation]]
and [[ Interpolation schemes | Runtime-Parameters:-Interpolation##INT_TABLE]] = 8.

<a name="SPEC_INT_GHOST_BOUNDARY"></a>
* #### `SPEC_INT_GHOST_BOUNDARY` &ensp; (&#8805;1) &ensp; [4]
    * **Description:**
Ghost boundary size for spectral interpolation.
    * **Restriction:**
Only applicable when the enabling compilation option
[[ --spectral_interpolation | Installation:-Option-List#--spectral_interpolation]]
and [[ Interpolation schemes | Runtime-Parameters:-Interpolation##INT_TABLE]] = 8.

## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of ELBDM | ELBDM]]