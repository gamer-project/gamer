Parameters described on this page:
[OPT__BC_FLU_XM](#OPT__BC_FLU_XM), &nbsp;
[OPT__BC_FLU_XP](#OPT__BC_FLU_XP), &nbsp;
[OPT__BC_FLU_YM](#OPT__BC_FLU_YM), &nbsp;
[OPT__BC_FLU_YP](#OPT__BC_FLU_YP), &nbsp;
[OPT__BC_FLU_ZM](#OPT__BC_FLU_ZM), &nbsp;
[OPT__BC_FLU_ZP](#OPT__BC_FLU_ZP), &nbsp;
[GAMMA](#GAMMA), &nbsp;
[MOLECULAR_WEIGHT](#MOLECULAR_WEIGHT), &nbsp;
[MU_NORM](#MU_NORM), &nbsp;
[ISO_TEMP](#ISO_TEMP), &nbsp;
[OPT__LR_LIMITER](#OPT__LR_LIMITER), &nbsp;
[MINMOD_COEFF](#MINMOD_COEFF), &nbsp;
[MINMOD_MAX_ITER](#MINMOD_MAX_ITER), &nbsp;
[OPT__1ST_FLUX_CORR](#OPT__1ST_FLUX_CORR), &nbsp;
[OPT__1ST_FLUX_CORR_SCHEME](#OPT__1ST_FLUX_CORR_SCHEME), &nbsp;
[DUAL_ENERGY_SWITCH](#DUAL_ENERGY_SWITCH), &nbsp;
[OPT__SAME_INTERFACE_B](#OPT__SAME_INTERFACE_B), &nbsp;
[OPT__FIXUP_FLUX](#OPT__FIXUP_FLUX), &nbsp;
[OPT__FIXUP_ELECTRIC](#OPT__FIXUP_ELECTRIC), &nbsp;
[OPT__FIXUP_RESTRICT](#OPT__FIXUP_RESTRICT), &nbsp;
[OPT__CORR_AFTER_ALL_SYNC](#OPT__CORR_AFTER_ALL_SYNC), &nbsp;
[OPT__NORMALIZE_PASSIVE](#OPT__NORMALIZE_PASSIVE), &nbsp;
[OPT__INT_FRAC_PASSIVE_LR](#OPT__INT_FRAC_PASSIVE_LR), &nbsp;
[OPT__RESET_FLUID](#OPT__RESET_FLUID), &nbsp;
[OPT__RESET_FLUID_INIT](#OPT__RESET_FLUID_INIT), &nbsp;
[OPT__FREEZE_FLUID](#OPT__FREEZE_FLUID), &nbsp;
[MIN_DENS](#MIN_DENS), &nbsp;
[MIN_PRES](#MIN_PRES), &nbsp;
[MIN_EINT](#MIN_EINT), &nbsp;
[OPT__CHECK_PRES_AFTER_FLU](#OPT__CHECK_PRES_AFTER_FLU), &nbsp;
[OPT__LAST_RESORT_FLOOR](#OPT__LAST_RESORT_FLOOR), &nbsp;
[JEANS_MIN_PRES](#JEANS_MIN_PRES), &nbsp;
[JEANS_MIN_PRES_LEVEL](#JEANS_MIN_PRES_LEVEL), &nbsp;
[JEANS_MIN_PRES_NCELL](#JEANS_MIN_PRES_NCELL), &nbsp;
[GAMMA_CR](#GAMMA_CR), &nbsp;
[CR_DIFF_PARA](#CR_DIFF_PARA), &nbsp;
[CR_DIFF_PERP](#CR_DIFF_PERP), &nbsp;
[CR_DIFF_MIN_B](#CR_DIFF_MIN_B) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="OPT__BC_FLU_XM"></a>
* #### `OPT__BC_FLU_XM` &ensp; (1=periodic, 2=outflow, 3=reflecting, 4=user, 5=diode) &ensp; [none]
    * **Description:**
Fluid boundary conditions on the -x face.
For the user-specified (i.e., inflow) boundary conditions, edit
`src/Fluid/Flu_BoundaryCondition_User.cpp`
for the cell-centered fluid variables and
`src/Model_Hydro/MHD_BoundaryCondition_User.cpp`
for the face-centered magnetic field, respectively,
when [[TESTPROB_ID | Runtime Parameters:-General#TESTPROB_ID]]=0
or a problem-specific function
when [[TESTPROB_ID | Runtime Parameters:-General#TESTPROB_ID]]&#8800;0.
    * **Restriction:**
When adopting periodic conditions, they must be applied to both the -x and +x faces.
Particles only support the periodic boundary conditions. When adopting any non-periodic
boundary condition for fluid, particles will be removed when getting too close
to the boundaries (see [[ PAR_REMOVE_CELL | Runtime-Parameters:-Particles#PAR_REMOVE_CELL ]] ).

<a name="OPT__BC_FLU_XP"></a>
* #### `OPT__BC_FLU_XP` &ensp; (1=periodic, 2=outflow, 3=reflecting, 4=user, 5=diode) &ensp; [none]
    * **Description:**
Fluid boundary conditions on the +x face.
See the description of [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).
    * **Restriction:**
See the restriction on [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).

<a name="OPT__BC_FLU_YM"></a>
* #### `OPT__BC_FLU_YM` &ensp; (1=periodic, 2=outflow, 3=reflecting, 4=user, 5=diode) &ensp; [none]
    * **Description:**
Fluid boundary conditions on the -y face.
See the description of [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).
    * **Restriction:**
See the restriction on [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).

<a name="OPT__BC_FLU_YP"></a>
* #### `OPT__BC_FLU_YP` &ensp; (1=periodic, 2=outflow, 3=reflecting, 4=user, 5=diode) &ensp; [none]
    * **Description:**
Fluid boundary conditions on the +y face.
See the description of [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).
    * **Restriction:**
See the restriction on [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).

<a name="OPT__BC_FLU_ZM"></a>
* #### `OPT__BC_FLU_ZM` &ensp; (1=periodic, 2=outflow, 3=reflecting, 4=user, 5=diode) &ensp; [none]
    * **Description:**
Fluid boundary conditions on the -z face.
See the description of [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).
    * **Restriction:**
See the restriction on [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).

<a name="OPT__BC_FLU_ZP"></a>
* #### `OPT__BC_FLU_ZP` &ensp; (1=periodic, 2=outflow, 3=reflecting, 4=user, 5=diode) &ensp; [none]
    * **Description:**
Fluid boundary conditions on the +z face.
See the description of [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).
    * **Restriction:**
See the restriction on [OPT__BC_FLU_XM](#OPT__BC_FLU_XM).

<a name="GAMMA"></a>
* #### `GAMMA` &ensp; (>1.0) &ensp; [5.0/3.0]
    * **Description:**
Ratio of specific heats (i.e., adiabatic index) in the ideal gas EoS.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--eos | Installation:-Option-List#--eos]]=`GAMMA`.
Does not support the isothermal EoS (i.e., `GAMMA=1.0`), for which
one should adopt [[--eos | Installation:-Option-List#--eos]]=`ISOTHERMAL`.

<a name="MOLECULAR_WEIGHT"></a>
* #### `MOLECULAR_WEIGHT` &ensp; (>0.0) &ensp; [0.6]
    * **Description:**
Mean molecular weight.
    * **Restriction:**

<a name="MU_NORM"></a>
* #### `MU_NORM` &ensp; (<0.0=m_H, 0.0=amu, >0.0=input manually in gram) &ensp; [-1.0]
    * **Description:**
Normalization of [MOLECULAR_WEIGHT](#MOLECULAR_WEIGHT).
    * **Restriction:**

<a name="ISO_TEMP"></a>
* #### `ISO_TEMP` &ensp; (>0.0) &ensp; [none]
    * **Description:**
Temperature in kelvin for the isothermal equation of state.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--eos | Installation:-Option-List#--eos]]=`ISOTHERMAL`.

<a name="OPT__LR_LIMITER"></a>
* #### `OPT__LR_LIMITER` &ensp; (-1&#8594; set to default, 0=none, 1=van Leer, 2=generalized minmod, 3=van Albada, 4=van Leer+generalized minmod, 6=central, 7=Athena) &ensp; [-1]
    * **Description:**
Slope limiter for data reconstruction. The coefficient of the
generalized minmod limiter can be set by [MINMOD_COEFF](#MINMOD_COEFF).
`7=Athena` mimics the extrema-preserving limiter implemented in Athena++.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--flu_scheme | Installation:-Option-List#--flu_scheme]]=`MHM/MHM_RP/CTU`.
`7=Athena` must work with
[[--slope | Installation:-Option-List#--slope]]=`PPM`.

<a name="MINMOD_COEFF"></a>
* #### `MINMOD_COEFF` &ensp; (1.0 &#8804; input &#8804; 2.0) &ensp; [1.5]
    * **Description:**
Coefficient of the generalized minmod limiter for data reconstruction
([OPT__LR_LIMITER](#OPT__LR_LIMITER)=2/4).
Smaller values lead to more diffusive but also more stable
solutions.
    * **Restriction:**

<a name="MINMOD_MAX_ITER"></a>
* #### `MINMOD_MAX_ITER` &ensp; (0=off, &#8805;0=on) &ensp; [0]
    * **Description:**
Maximum number of iterations to reduce [MINMOD_COEFF](#MINMOD_COEFF)
when data reconstruction fails. It improves code stability but
may break strict conservation.
    * **Restriction:**

<a name="OPT__1ST_FLUX_CORR"></a>
* #### `OPT__1ST_FLUX_CORR` &ensp; (<0 &#8594; set to default, 0=off, 1=3D, 2=3D+1D) &ensp; [2]
    * **Description:**
Correct unphysical results (e.g., negative density or pressure) by
recalculating the solutions of failed cells with only first-order
accuracy in space and time.
`OPT__1ST_FLUX_CORR=1`: only use dimensionally unsplit update.
`OPT__1ST_FLUX_CORR=2`: use dimensionally unsplit update first;
and if it fails, switch to dimensionally split update.
The number of cells corrected by this option will be recorded in the file
[[Record__NCorrUnphy | Simulation-Logs:-Record__NCorrUnphy]].
    * **Restriction:**
Only applicable when adopting the compilation option
[[--flu_scheme | Installation:-Option-List#--flu_scheme]]=`MHM/MHM_RP/CTU`.
[[--mhd | Installation:-Option-List#--mhd]] currently does not support `3D+1D`.
Be aware that this option may cause conservation errors.

<a name="OPT__1ST_FLUX_CORR_SCHEME"></a>
* #### `OPT__1ST_FLUX_CORR_SCHEME` &ensp; (<0 &#8594; set to default, 0=none, 1=Roe, 2=HLLC, 3=HLLE, 4=HLLD) &ensp; [1]
    * **Description:**
Riemann solver for `OPT__1ST_FLUX_CORR`.
    * **Restriction:**

<a name="DUAL_ENERGY_SWITCH"></a>
* #### `DUAL_ENERGY_SWITCH` &ensp; (&#8805;0.0) &ensp; [2e-2]
    * **Description:**
Switch threshold of the dual energy formalism. Specifically, we enable
the dual energy formalism for cells with
<var>E</var><sub><var>int<var><var></sub>/<var>E</var><sub><var>kin<var><var></sub>
&#8287;&lt;&#8287; <var>&xi;</var>, where
<var>E</var><sub><var>int<var></sub> is gas internal energy,
<var>E</var><sub><var>kin<var></sub> is gas kinematic energy,
and <var>&xi;</var>=`DUAL_ENERGY_SWITCH`.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--dual | Installation:-Option-List#--dual]].

<a name="OPT__SAME_INTERFACE_B"></a>
* #### `OPT__SAME_INTERFACE_B` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Ensure adjacent patches on the same level have exactly the same magnetic field
on their shared interfaces (including round-off errors).
This option is mainly for debugging purposes.
    * **Restriction:**

<a name="OPT__FIXUP_FLUX"></a>
* #### `OPT__FIXUP_FLUX` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Correct the coarse-grid data by the fluxes on the coarse-fine boundaries.
    * **Restriction:**

<a name="OPT__FIXUP_ELECTRIC"></a>
* #### `OPT__FIXUP_ELECTRIC` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Correct the coarse-grid magnetic field by the electric field on the
coarse-fine boundaries.
    * **Restriction:**
For [[--mhd | Installation:-Option-List#--mhd]] only.

<a name="OPT__FIXUP_RESTRICT"></a>
* #### `OPT__FIXUP_RESTRICT` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Correct the coarse-grid data by the volume-weighted average of
the child patch data.
    * **Restriction:**

<a name="OPT__CORR_AFTER_ALL_SYNC"></a>
* #### `OPT__CORR_AFTER_ALL_SYNC` &ensp; (-1=set to default, 0=off, 1=every step, 2=before dump) &ensp; [depend]
    * **Description:**
Apply additional corrections when all levels are synchronized.
See `Fluid/Flu_CorrAfterAllSync.cpp` for details. This functionality
is mainly for achieving [[bitwise reproducibility | Bitwise Reproducibility]].
`OPT__CORR_AFTER_ALL_SYNC=1`: apply corrections after each root-level update.
`OPT__CORR_AFTER_ALL_SYNC=2`: apply corrections before each data dump.
The default depends on the compilation option
[[--debug | Installation:-Option-List#--debug]].
    * **Restriction:**
Must be turned on when enabling the compilation option
[[--bitwise_reproducibility | Installation:-Option-List#--bitwise_reproducibility]].

<a name="OPT__NORMALIZE_PASSIVE"></a>
* #### `OPT__NORMALIZE_PASSIVE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Renormalize passively advected scalars after every update to ensure that
the sum of their mass fractions equals unity. See also
[[ Add Problem-specific Grid Fields and Particle Attributes | Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes ]].
    * **Restriction:**

<a name="OPT__INT_FRAC_PASSIVE_LR"></a>
* #### `OPT__INT_FRAC_PASSIVE_LR` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Convert specified passive scalars from mass density to fraction during data reconstruction. See also
[[ Add Problem-specific Grid Fields and Particle Attributes | Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes ]].
    * **Restriction:**

<a name="OPT__RESET_FLUID"></a>
* #### `OPT__RESET_FLUID` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Reset fluid variables after each update. It can be useful for adding
sources or sinks to a simulation. To specify how to reset fluid variables,
edit the source file `src/Fluid/Flu_ResetByUser.cpp` directly or add a
specific routine to your test problem files. See [[Adding New Simulations]]
for details.
    * **Restriction:**

<a name="OPT__RESET_FLUID_INIT"></a>
* #### `OPT__RESET_FLUID_INIT` &ensp; (<0 &#8594; set to [OPT__RESET_FLUID](#OPT__RESET_FLUID), 0=off, 1=on) &ensp; [-1]
    * **Description:**
Reset fluid variables during initialization. See [OPT__RESET_FLUID](#OPT__RESET_FLUID) for details.
    * **Restriction:**

<a name="OPT__FREEZE_FLUID"></a>
* #### `OPT__FREEZE_FLUID` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Do not evolve fluid variables at all. It can be useful for evolving particles
in a static gravitational potential of fluid or wave dark matter.
    * **Restriction:**

<a name="OPT__CHECK_PRES_AFTER_FLU"></a>
* #### `OPT__CHECK_PRES_AFTER_FLU` &ensp; (<0 &#8594; set to default, 0=off, 1=on) &ensp; [-1]
    * **Description:**
Check unphysical pressure at the end of the fluid solver. If it is off,
the code will only check internal energy. This is particularly useful
for a complicated equation of state.
    * **Restriction:**

<a name="OPT__LAST_RESORT_FLOOR"></a>
* #### `OPT__LAST_RESORT_FLOOR` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Apply floor values (e.g., [MIN_DENS](#MIN_DENS) and [MIN_EINT](#MIN_EINT)) when both
[OPT__1ST_FLUX_CORR](#OPT__1ST_FLUX_CORR) and
[[AUTO_REDUCE_DT | Runtime Parameters:-Timestep#AUTO_REDUCE_DT]]
fail. Always enable this option unless for debugging purposes.
    * **Restriction:**

<a name="MIN_DENS"></a>
* #### `MIN_DENS` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Gas mass density floor.
    * **Restriction:**

<a name="MIN_PRES"></a>
* #### `MIN_PRES` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Gas pressure floor.
    * **Restriction:**

<a name="MIN_EINT"></a>
* #### `MIN_EINT` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Gas internal energy density floor.
    * **Restriction:**

<a name="JEANS_MIN_PRES"></a>
* #### `JEANS_MIN_PRES` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Employ a pressure floor in the hydrodynamic solver to ensure that the
local Jeans length of all cells on all levels is resolved by at least
a width equal to
[JEANS_MIN_PRES_NCELL](#JEANS_MIN_PRES_NCELL) cells on level
[JEANS_MIN_PRES_LEVEL](#JEANS_MIN_PRES_LEVEL).
    * **Restriction:**
Only applicable when enabling the compilation option
[[--gravity | Installation:-Option-List#--gravity]].

<a name="JEANS_MIN_PRES_LEVEL"></a>
* #### `JEANS_MIN_PRES_LEVEL` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [ [[MAX_LEVEL | Runtime Parameters:-Refinement#MAX_LEVEL]] ]
    * **Description:**
See [JEANS_MIN_PRES](#JEANS_MIN_PRES).
    * **Restriction:**

<a name="JEANS_MIN_PRES_NCELL"></a>
* #### `JEANS_MIN_PRES_NCELL` &ensp; (&#8805;1) &ensp; [4]
    * **Description:**
See [JEANS_MIN_PRES](#JEANS_MIN_PRES).
    * **Restriction:**
Must be an integer.

<a name="GAMMA_CR"></a>
* #### `GAMMA_CR` &ensp; (>1.0) &ensp; [4.0/3.0]
    * **Description:**
Effective adiabatic index of cosmic rays.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--cosmic_ray | Installation:-Option-List#--cosmic_ray]].

<a name="CR_DIFF_PARA"></a>
* #### `CR_DIFF_PARA` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Cosmic-ray diffusion coefficient parallel to the magnetic field.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--cr_diffusion | Installation:-Option-List#--cr_diffusion]].

<a name="CR_DIFF_PERP"></a>
* #### `CR_DIFF_PERP` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Cosmic-ray diffusion coefficient perpendicular to the magnetic field.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--cr_diffusion | Installation:-Option-List#--cr_diffusion]].

<a name="CR_DIFF_MIN_B"></a>
* #### `CR_DIFF_MIN_B` &ensp; (none) &ensp; [0.0]
    * **Description:**
Disable cosmic-ray diffusion locally when the magnetic field amplitude is smaller
than this threshold.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--cr_diffusion | Installation:-Option-List#--cr_diffusion]].


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Hydro | Hydro]]
