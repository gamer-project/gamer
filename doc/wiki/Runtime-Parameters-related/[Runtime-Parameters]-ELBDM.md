
## Compilation Options

Related options:
[[--model=ELBDM | [Installation]-Option-List#--model]] &nbsp;

## Runtime Parameters

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
[ELBDM_RESCALE_MASS_ERROR](#ELBDM_RESCALE_MASS_ERROR), &nbsp;
[ELBDM_RESCALE_MASS_STEPS](#ELBDM_RESCALE_MASS_STEPS), &nbsp;
[OPT__INT_PHASE](#OPT__INT_PHASE), &nbsp;
[OPT__RES_PHASE](#OPT__RES_PHASE), &nbsp;
[SPEC_INT_TABLE_PATH](#SPEC_INT_TABLE_PATH), &nbsp;
[SPEC_INT_XY_INSTEAD_DEPHA](#SPEC_INT_XY_INSTEAD_DEPHA), &nbsp;
[SPEC_INT_VORTEX_THRESHOLD](#SPEC_INT_VORTEX_THRESHOLD), &nbsp;
[SPEC_INT_GHOST_BOUNDARY](#SPEC_INT_GHOST_BOUNDARY) &nbsp;

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="ELBDM_MASS"></a>
* #### `ELBDM_MASS` &ensp; (>0.0) &ensp; [none]
    * **Description:**
Particle mass in $eV/c^2$.
Note that the input unit is fixed regardless of whether
[[OPT__UNIT | [Runtime-Parameters]-Units#OPT__UNIT]]
or
[[--comoving | [Installation]-Option-List#--comoving]]
is enabled.
    * **Restriction:**

<a name="ELBDM_PLANCK_CONST"></a>
* #### `ELBDM_PLANCK_CONST` &ensp; (>0.0) &ensp; [conform to the unit system set by [[OPT__UNIT | [Runtime-Parameters]-Units#OPT__UNIT]] or [[--comoving | [Installation]-Option-List#--comoving]]]
    * **Description:**
Reduced Planck constant in $g\ cm^2/s$.
    * **Restriction:**
The input value will be overwritten by the default value when
[[OPT__UNIT | [Runtime-Parameters]-Units#OPT__UNIT]]
or
[[--comoving | [Installation]-Option-List#--comoving]]
is enabled.

<a name="ELBDM_LAMBDA"></a>
* #### `ELBDM_LAMBDA` &ensp; (any floating value) &ensp; [1.0]
    * **Description:**
Quartic self-interaction coefficient in ELBDM.
    * **Restriction:**
Only applicable when the compilation option
[[--self_interaction | [Installation]-Option-List#--self_interaction]] is enabled.

<a name="ELBDM_TAYLOR3_COEFF"></a>
* #### `ELBDM_TAYLOR3_COEFF` &ensp; (&#8805;0.125) &ensp; [1.0/6.0]
    * **Description:**
Coefficient of the 3rd-order term in the Taylor expansion for the finite-difference wave solver.
    * **Restriction:**
Only applicable when the compilation option
[[--wave_scheme | [Installation]-Option-List#--wave_scheme]]=`FD`
is enabled.
Ignored if [ELBDM_TAYLOR3_AUTO](#ELBDM_TAYLOR3_AUTO) is enabled.
Values &#8804; 0.125 are always unstable.
Values &#8804; 1/6 are unstable if
[[DT__FLUID | [Runtime-Parameters]-Timestep#DT__FLUID]] > $\sqrt{27}\pi/32$
(or $\sqrt{3}\pi/8$)
when [[--laplacian_four | [Installation]-Option-List#--laplacian_four]]
is enabled (or disabled), respectively.

<a name="ELBDM_TAYLOR3_AUTO"></a>
* #### `ELBDM_TAYLOR3_AUTO` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Automatically determine
[ELBDM_TAYLOR3_COEFF](#ELBDM_TAYLOR3_COEFF) to minimize the amplitude error
at the smallest wavelength.
    * **Restriction:**
Useless if [[OPT__FREEZE_FLUID | [Runtime-Parameters]-Hydro#OPT__FREEZE_FLUID]] is on.

<a name="ELBDM_REMOVE_MOTION_CM"></a>
* #### `ELBDM_REMOVE_MOTION_CM` &ensp; (0=none, 1=init, 2=every step) &ensp; [0]
    * **Description:**
Remove the center-of-mass velocity.
    * **Restriction:**
Only applicable when
[[OPT__CK_CONSERVATION | [Runtime-Parameters]-Miscellaneous#OPT__CK_CONSERVATION]]
is enabled.
Not supported when
[[--bitwise_reproducibility | [Installation]-Option-List#--bitwise_reproducibility]]=true.

<a name="ELBDM_BASE_SPECTRAL"></a>
* #### `ELBDM_BASE_SPECTRAL` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Adopt the spectral method to evolve the base-level wave function.
    * **Restriction:**
Requires [[--fftw | [Installation]-Option-List#--fftw]]=FFTW2/FFTW3
and periodic boundary conditions in all directions:
[[OPT__BC_FLU | [Runtime-Parameters]-Hydro#OPT__BC_FLU_XM]]=1.

<a name="ELBDM_MATCH_PHASE"></a>
* #### `ELBDM_MATCH_PHASE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
During data restriction, unwrap the average phases of child patches on a wave level
to match the phases of their corresponding parent patches on a fluid level.
    * **Restriction:**
Only applicable when enabling the compilation option
[[ --elbdm_scheme | [Installation]-Option-List#--elbdm_scheme]]=`HYBRID`.
Requires [[ OPT__UM_IC_LEVEL | [Runtime-Parameters]-Initial-Conditions#OPT__UM_IC_LEVEL ]]
< [ELBDM_FIRST_WAVE_LEVEL](#ELBDM_FIRST_WAVE_LEVEL).

<a name="ELBDM_FIRST_WAVE_LEVEL"></a>
* #### `ELBDM_FIRST_WAVE_LEVEL` &ensp; (1 &#8804; input &#8804; [[ MAX_LEVEL | [Runtime-Parameters]-Refinement#MAX_LEVEL]]) &ensp; [none]
    * **Description:**
Level at which to switch to the wave solver.
    * **Restriction:**
Only applicable when enabling the compilation option
[[ --elbdm_scheme | [Installation]-Option-List#--elbdm_scheme]]=`HYBRID`.

<a name="ELBDM_RESCALE_MASS_ERROR"></a>
* #### `ELBDM_RESCALE_MASS_ERROR` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Rescale the total ELBDM mass to its initial value every [ELBDM_RESCALE_MASS_STEPS](#ELBDM_RESCALE_MASS_STEPS) steps
to ensure mass conservation.
    * **Restriction:**
Only applicable when enabling
[[ OPT__CK_CONSERVATION | [Runtime-Parameters]-Miscellaneous#OPT__CK_CONSERVATION ]].

<a name="ELBDM_RESCALE_MASS_STEPS"></a>
* #### `ELBDM_RESCALE_MASS_STEPS` &ensp; (&#8805;1) &ensp; [100]
    * **Description:**
See [ELBDM_RESCALE_MASS_ERROR](#ELBDM_RESCALE_MASS_ERROR).
    * **Restriction:**
Only applicable when enabling [ELBDM_RESCALE_MASS_ERROR](#ELBDM_RESCALE_MASS_ERROR).

<a name="OPT__INT_PHASE"></a>
* #### `OPT__INT_PHASE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Perform data interpolation on the phase field rather than the wave function itself
when both the parent and child patches are on wave levels.
The interpolation scheme is determined by
[[OPT__FLU_INT_SCHEME | [Runtime-Parameters]-Interpolation#OPT__FLU_INT_SCHEME]] and
[[OPT__REF_FLU_INT_SCHEME | [Runtime-Parameters]-Interpolation#OPT__REF_FLU_INT_SCHEME]].
See also [OPT__RES_PHASE](#OPT__RES_PHASE).
    * **Restriction:**
The "1D MinMod limiter" interpolation scheme is not supported.

<a name="OPT__RES_PHASE"></a>
* #### `OPT__RES_PHASE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Perform data restriction on the phase field rather than the wave function itself
when both the parent and child patches are on wave levels.
See also [OPT__INT_PHASE](#OPT__INT_PHASE).
    * **Restriction:**

<a name="SPEC_INT_TABLE_PATH"></a>
* #### `SPEC_INT_TABLE_PATH` &ensp; (none) &ensp; [none]
    * **Description:**
Path to the spectral interpolation table.
See [[ELBDM Spectral Interpolation | [ELBDM]-Spectral-Interpolation]] for details.
A script for downloading the table is available at
`example/test_problem/ELBDM/LSS_Hybrid/download_spectral_interpolation_tables.sh`.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--spectral_interpolation | [Installation]-Option-List#--spectral_interpolation]]
and adopting [[Interpolation Scheme | [Runtime-Parameters]-Interpolation]]=8.

<a name="SPEC_INT_XY_INSTEAD_DEPHA"></a>
* #### `SPEC_INT_XY_INSTEAD_DEPHA` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Interpolate x and y (real and imaginary parts in current implementation) around vortices
instead of density and phase for the spectral interpolation,
which has the advantage of being well-defined across vortices.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--spectral_interpolation | [Installation]-Option-List#--spectral_interpolation]]
and adopting [[Interpolation Scheme | [Runtime-Parameters]-Interpolation]]=8.

<a name="SPEC_INT_VORTEX_THRESHOLD"></a>
* #### `SPEC_INT_VORTEX_THRESHOLD` &ensp; (&#8805;0.0) &ensp; [0.1]
    * **Description:**
Vortex detection threshold for [SPEC_INT_XY_INSTEAD_DEPHA](#SPEC_INT_XY_INSTEAD_DEPHA),
triggered when $\nabla^2 S\ dx^2 > \rm threshold$, indicating a significant phase jump.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--spectral_interpolation | [Installation]-Option-List#--spectral_interpolation]]
and adopting [[Interpolation Scheme | [Runtime-Parameters]-Interpolation]]=8.

<a name="SPEC_INT_GHOST_BOUNDARY"></a>
* #### `SPEC_INT_GHOST_BOUNDARY` &ensp; (&#8805;1) &ensp; [4]
    * **Description:**
Ghost boundary size for the spectral interpolation.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--spectral_interpolation | [Installation]-Option-List#--spectral_interpolation]]
and adopting [[Interpolation Scheme | [Runtime-Parameters]-Interpolation]]=8.

## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of ELBDM | ELBDM]]