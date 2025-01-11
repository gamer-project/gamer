## Enabling AMR

It only takes three steps to enable AMR:

* Set [MAX_LEVEL](#MAX_LEVEL)
* Turn on at least one of the refinement criteria `OPT__FLAG_*`
* Edit the corresponding input file(s)
[[Input__Flag_* | Runtime-Parameters:-Input__Flag_*]]
to specify the refinement thresholds

See the descriptions of various refinement criteria `OPT__FLAG_*`
given on this page for details.


## Compilation Options

Related options:
[[--nlevel | Installation:-Option-List#--nlevel]], &nbsp;
[[--max_patch | Installation:-Option-List#--max_patch]] &nbsp;


## Runtime Parameters

Parameters described on this page:
[REGRID_COUNT](#REGRID_COUNT), &nbsp;
[REFINE_NLEVEL](#REFINE_NLEVEL), &nbsp;
[FLAG_BUFFER_SIZE](#FLAG_BUFFER_SIZE), &nbsp;
[FLAG_BUFFER_SIZE_MAXM1_LV](#FLAG_BUFFER_SIZE_MAXM1_LV), &nbsp;
[FLAG_BUFFER_SIZE_MAXM2_LV](#FLAG_BUFFER_SIZE_MAXM2_LV), &nbsp;
[MAX_LEVEL](#MAX_LEVEL), &nbsp;
[OPT__FLAG_RHO](#OPT__FLAG_RHO), &nbsp;
[OPT__FLAG_RHO_GRADIENT](#OPT__FLAG_RHO_GRADIENT), &nbsp;
[OPT__FLAG_PRES_GRADIENT](#OPT__FLAG_PRES_GRADIENT), &nbsp;
[OPT__FLAG_LRTZ_GRADIENT](#OPT__FLAG_LRTZ_GRADIENT), &nbsp;
[OPT__FLAG_VORTICITY](#OPT__FLAG_VORTICITY), &nbsp;
[OPT__FLAG_JEANS](#OPT__FLAG_JEANS), &nbsp;
[OPT__FLAG_CURRENT](#OPT__FLAG_CURRENT), &nbsp;
[OPT__FLAG_CRAY](#OPT__FLAG_CRAY), &nbsp;
[OPT__FLAG_LOHNER_DENS](#OPT__FLAG_LOHNER_DENS), &nbsp;
[OPT__FLAG_LOHNER_ENGY](#OPT__FLAG_LOHNER_ENGY), &nbsp;
[OPT__FLAG_LOHNER_PRES](#OPT__FLAG_LOHNER_PRES), &nbsp;
[OPT__FLAG_LOHNER_TEMP](#OPT__FLAG_LOHNER_TEMP), &nbsp;
[OPT__FLAG_LOHNER_ENTR](#OPT__FLAG_LOHNER_ENTR), &nbsp;
[OPT__FLAG_LOHNER_CRAY](#OPT__FLAG_LOHNER_CRAY), &nbsp;
[OPT__FLAG_LOHNER_FORM](#OPT__FLAG_LOHNER_FORM), &nbsp;
[OPT__FLAG_USER](#OPT__FLAG_USER), &nbsp;
[OPT__FLAG_USER_NUM](#OPT__FLAG_USER_NUM), &nbsp;
[OPT__FLAG_REGION](#OPT__FLAG_REGION), &nbsp;
[OPT__FLAG_ANGULAR](#OPT__FLAG_ANGULAR), &nbsp;
[FLAG_ANGULAR_CEN_X](#FLAG_ANGULAR_CEN_X), &nbsp;
[FLAG_ANGULAR_CEN_Y](#FLAG_ANGULAR_CEN_Y), &nbsp;
[FLAG_ANGULAR_CEN_Z](#FLAG_ANGULAR_CEN_Z), &nbsp;
[OPT__FLAG_RADIAL](#OPT__FLAG_RADIAL), &nbsp;
[FLAG_RADIAL_CEN_X](#FLAG_RADIAL_CEN_X), &nbsp;
[FLAG_RADIAL_CEN_Y](#FLAG_RADIAL_CEN_Y), &nbsp;
[FLAG_RADIAL_CEN_Z](#FLAG_RADIAL_CEN_Z), &nbsp;
[OPT__FLAG_NPAR_PATCH](#OPT__FLAG_NPAR_PATCH), &nbsp;
[OPT__FLAG_NPAR_CELL](#OPT__FLAG_NPAR_CELL), &nbsp;
[OPT__FLAG_PAR_MASS_CELL](#OPT__FLAG_PAR_MASS_CELL), &nbsp;
[OPT__NO_FLAG_NEAR_BOUNDARY](#OPT__NO_FLAG_NEAR_BOUNDARY), &nbsp;
[OPT__PATCH_COUNT](#OPT__PATCH_COUNT), &nbsp;
[OPT__PARTICLE_COUNT](#OPT__PARTICLE_COUNT), &nbsp;
[OPT__REUSE_MEMORY](#OPT__REUSE_MEMORY), &nbsp;
[OPT__MEMORY_POOL](#OPT__MEMORY_POOL) &nbsp;

Other related parameters:
[[OPT__UM_IC_DOWNGRADE | Runtime-Parameters:-Initial-Conditions#OPT__UM_IC_DOWNGRADE]], &nbsp;
[[OPT__UM_IC_REFINE | Runtime-Parameters:-Initial-Conditions#OPT__UM_IC_REFINE]], &nbsp;
[[OPT__CK_REFINE | Runtime Parameters:-Miscellaneous#OPT__CK_REFINE]], &nbsp;
[[Interpolation Schemes | Runtime Parameters:-Interpolation]] &nbsp;

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="REGRID_COUNT"></a>
* #### `REGRID_COUNT` &ensp; (&#8805;1) &ensp; [4]
    * **Description:**
Check refinement on level <var>l</var> to create and remove patches on
level <var>l+1</var> every `REGRID_COUNT` substeps on level <var>l</var>.
Adopting a larger value can improve performance but may also increase the
risk of allowing the phenomenon of interest to leave regions with adequate
resolution. Setting this parameter to an extremely large value is essentially
the same as conducting a static mesh refinement (SMR) simulation.
    * **Restriction:**

<a name="REFINE_NLEVEL"></a>
* #### `REFINE_NLEVEL` &ensp; (&#8805;1) &ensp; [1]
    * **Description:**
Number of new AMR levels to be created at once during refinement.
    * **Restriction:**

<a name="FLAG_BUFFER_SIZE"></a>
* #### `FLAG_BUFFER_SIZE` &ensp; (0 &#8804; input &#8804; [[--patch_size | Installation:-Option-List#--patch_size]]) &ensp; [ [[--patch_size | Installation:-Option-List#--patch_size]] ]
    * **Description:**
Number of flag buffers (denoted as <var>N</var><sub>buf</sub>).
When checking refinement, a patch is flagged for refinement if any of
its cells satisfy the refinement criteria. In addition, we add
<var>(1+2N</var><sub>buf</sub><var>)</var><sup>3</sup><var>-1</var>
flag buffers around each flagged cell and refine the corresponding
sibling patches as well if any of these flag buffers extend across
the patch border. See
[Fig. 2 in Schive et al. 2010](http://iopscience.iop.org/article/10.1088/0067-0049/186/2/457/meta#apjs325434f2)
for an illustration.
Adopting a value smaller than `PATCH_SIZE` may increase the risk of allowing the
phenomenon of interest leaving regions with adequate resolution.

        Note that to alleviate the issue of over-refinemenet resulting
from a large <var>N</var><sub>buf</sub> (e.g., `PATCH_SIZE`), we only apply this
parameter to levels below [MAX_LEVEL](#MAX_LEVEL)-2.
<var>N</var><sub>buf</sub> on levels [MAX_LEVEL](#MAX_LEVEL)-1 and
[MAX_LEVEL](#MAX_LEVEL)-2 are set by
[FLAG_BUFFER_SIZE_MAXM1_LV](#FLAG_BUFFER_SIZE_MAXM1_LV) and
[FLAG_BUFFER_SIZE_MAXM2_LV](#FLAG_BUFFER_SIZE_MAXM2_LV), respectively,
which can be safely set to smaller values.
    * **Restriction:**

<a name="FLAG_BUFFER_SIZE_MAXM1_LV"></a>
* #### `FLAG_BUFFER_SIZE_MAXM1_LV` &ensp; (0 &#8804; input &#8804; [[--patch_size | Installation:-Option-List#--patch_size]]; <0 &#8594; set to default) &ensp; [ [REGRID_COUNT](#REGRID_COUNT) ]
    * **Description:**
<var>N</var><sub>buf</sub> on level [MAX_LEVEL](#MAX_LEVEL)-1. See
[FLAG_BUFFER_SIZE](#FLAG_BUFFER_SIZE) for details.
It is recommended to set this parameter to 2 ~ 4 (in accordance with [REGRID_COUNT](#REGRID_COUNT))
for better performance.
    * **Restriction:**

<a name="FLAG_BUFFER_SIZE_MAXM2_LV"></a>
* #### `FLAG_BUFFER_SIZE_MAXM2_LV` &ensp; (0 &#8804; input &#8804; [[--patch_size | Installation:-Option-List#--patch_size]]; <0 &#8594; set to default) &ensp; [depend]
    * **Description:**
<var>N</var><sub>buf</sub> on level [MAX_LEVEL](#MAX_LEVEL)-2. See
[FLAG_BUFFER_SIZE](#FLAG_BUFFER_SIZE) for details. The default value
depends on both [[--patch_size | Installation:-Option-List#--patch_size]] and
[FLAG_BUFFER_SIZE_MAXM1_LV](#FLAG_BUFFER_SIZE_MAXM1_LV).
    * **Restriction:**

<a name="MAX_LEVEL"></a>
* #### `MAX_LEVEL` &ensp; (0 &#8804; input < [[--nlevel | Installation:-Option-List#--nlevel]]) &ensp; [ [[--nlevel | Installation:-Option-List#--nlevel]]-1 ]
    * **Description:**
Maximum allowed AMR level. Do not confuse with the compilation option
[[--nlevel | Installation:-Option-List#--nlevel]]. One can regard
`MAX_LEVEL` as a runtime refinement criterion that simply forbids
creating any patch above that level, and `--nlevel - 1` is the upper
bound of `MAX_LEVEL`.
    * **Restriction:**

<a name="OPT__FLAG_RHO"></a>
* #### `OPT__FLAG_RHO` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: gas mass density. Specify the refinement thresholds
on different levels in the input file `Input__Flag_Rho` with the
[[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_Rho8`
(must rename it as `Input__Flag_Rho` to actually use it).
By setting the density threshold ratio between adjacent
levels to 8, it leads to a roughly fixed mass resolution similar to a
Lagrangian method.
    * **Restriction:**

<a name="OPT__FLAG_RHO_GRADIENT"></a>
* #### `OPT__FLAG_RHO_GRADIENT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: gas mass density gradient. Specifically, a cell
on level <var>l</var> will be flagged for refinement if its normalized
density slope along a spatial direction &xi; (= x/y/z
in Cartesian coordinates) satisfies
&#8706;<sub>&xi;</sub>&rho;&#8901;&Delta;&xi;<sub>l</sub>&#8287;/&#8287;&rho;&#8287;&#8805;&#8287;&eta;<sub>l</sub>,
where &rho; is gas mass density, &Delta;&xi;<sub>l</sub> is the cell
width along &xi; on level <var>l</var>, and &eta;<sub>l</sub> is the
refinement threshold on level <var>l</var>.
Specify the refinement thresholds on different levels in the input
file `Input__Flag_RhoGradient` with the
[[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_RhoGradient`.
    * **Restriction:**

<a name="OPT__FLAG_PRES_GRADIENT"></a>
* #### `OPT__FLAG_PRES_GRADIENT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: gas pressure gradient. See
[OPT__FLAG_RHO_GRADIENT](#OPT__FLAG_RHO_GRADIENT) for the definition
of normalized slope. Specify the refinement thresholds on different
levels in the input file `Input__Flag_PresGradient` with the
[[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_PresGradient`.
    * **Restriction:**

<a name="OPT__FLAG_LRTZ_GRADIENT"></a>
* #### `OPT__FLAG_LRTZ_GRADIENT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: Lorentz factor gradient. See
[OPT__FLAG_RHO_GRADIENT](#OPT__FLAG_RHO_GRADIENT) for the definition
of normalized slope. Specify the refinement thresholds on different
levels in the input file `Input__Flag_LrtzGradient` with the
[[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_LrtzGradient`.
    * **Restriction:**
Must compile with [[--srhd | Installation:-Option-List#--srhd]].

<a name="OPT__FLAG_VORTICITY"></a>
* #### `OPT__FLAG_VORTICITY` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: gas vorticity. Specifically, a cell
on level <var>l</var> will be flagged for refinement if its normalized
vorticity satisfies
|&#8711;&#10799;<var>v</var>|&#8901;&Delta;&xi;<sub>l</sub>&#8287;/&#8287;|<var>v</var>|&#8287;&#8805;&#8287;&eta;<sub>l</sub>,
where <var>v</var> is gas velocity, &Delta;&xi;<sub>l</sub> is the cell
width along &xi; on level <var>l</var>, and &eta;<sub>l</sub> is the
refinement threshold on level <var>l</var>. Specify the refinement
thresholds on different levels in the input file `Input__Flag_Vorticity`
with the [[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_Vorticity`.
    * **Restriction:**

<a name="OPT__FLAG_JEANS"></a>
* #### `OPT__FLAG_JEANS` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: gas Jeans length. It ensures that the Jeans length
is resolved by at least <var>N</var> cells. Specifically, a cell
on level <var>l</var> will be flagged for refinement if its estimated
Jeans length &lambda;<sub>J</sub> satisfies
&lambda;<sub>J</sub>&#8287;&#8801;&#8287;(&pi;&gamma;<var>P</var>/<var>G</var>&rho;<sup>2</sup>)<sup>1/2</sup>&#8287;<&#8287;<var>N</var><sub>l</sub>&Delta;&xi;<sub>l</sub>,
where &gamma; is adiabatic index ([[GAMMA | Runtime-Parameters:-Hydro#GAMMA]]), <var>P</var> is gas pressure, &rho; is gas mass density, <var>G</var> is gravitational constant, &Delta;&xi;<sub>l</sub> is the cell width along &xi; on level <var>l</var>, and <var>N</var><sub>l</sub> is the refinement threshold on level <var>l</var>.
Specify the refinement
thresholds on different levels in the input file `Input__Flag_Jeans`
with the [[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_Jeans`.
Recommended values: &#8805;4.
    * **Restriction:**

<a name="OPT__FLAG_CURRENT"></a>
* #### `OPT__FLAG_CURRENT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: current density in MHD. Specifically, a cell
on level <var>l</var> will be flagged for refinement if its currenty
density satisfies
|&#8711;&#10799;<var>B</var>|&#8901;&Delta;&xi;<sub>l</sub>&#8287;/&#8287;|<var>B</var>|&#8287;&#8805;&#8287;&eta;<sub>l</sub>,
where <var>B</var> is the magnetic field, &Delta;&xi;<sub>l</sub> is the cell
width along &xi; on level <var>l</var>, and &eta;<sub>l</sub> is the
refinement threshold on level <var>l</var>. Specify the refinement
thresholds on different levels in the input file `Input__Flag_Current`
with the [[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_Current`.
    * **Restriction:**
Must compile with [[--mhd | Installation:-Option-List#--mhd]].

<a name="OPT__FLAG_CRAY"></a>
* #### `OPT__FLAG_CRAY` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: cosmic-ray energy density. Specify the refinement
thresholds on different levels in the input file `Input__Flag_CRay`
with the [[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_CRay`.
    * **Restriction:**
Must compile with [[--cosmic_ray | Installation:-Option-List#--cosmic_ray]].

<a name="OPT__FLAG_LOHNER_DENS"></a>
* #### `OPT__FLAG_LOHNER_DENS` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: Loehner's error estimator on gas mass density.
See [OPT__FLAG_LOHNER_FORM](#OPT__FLAG_LOHNER_FORM) for details.
    * **Restriction:**

<a name="OPT__FLAG_LOHNER_ENGY"></a>
* #### `OPT__FLAG_LOHNER_ENGY` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: Loehner's error estimator on gas energy density
(excluding potential energy).
See [OPT__FLAG_LOHNER_FORM](#OPT__FLAG_LOHNER_FORM) for details.
    * **Restriction:**

<a name="OPT__FLAG_LOHNER_PRES"></a>
* #### `OPT__FLAG_LOHNER_PRES` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: Loehner's error estimator on gas pressure.
See [OPT__FLAG_LOHNER_FORM](#OPT__FLAG_LOHNER_FORM) for details.
    * **Restriction:**

<a name="OPT__FLAG_LOHNER_TEMP"></a>
* #### `OPT__FLAG_LOHNER_TEMP` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: Loehner's error estimator on gas temperature.
See [OPT__FLAG_LOHNER_FORM](#OPT__FLAG_LOHNER_FORM) for details.
    * **Restriction:**

<a name="OPT__FLAG_LOHNER_ENTR"></a>
* #### `OPT__FLAG_LOHNER_ENTR` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: Loehner's error estimator on gas entropy.
See [OPT__FLAG_LOHNER_FORM](#OPT__FLAG_LOHNER_FORM) for details.
    * **Restriction:**

<a name="OPT__FLAG_LOHNER_CRAY"></a>
* #### `OPT__FLAG_LOHNER_CRAY` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: Loehner's error estimator on cosmic-ray energy density.
See [OPT__FLAG_LOHNER_FORM](#OPT__FLAG_LOHNER_FORM) for details.
    * **Restriction:**
Must compile with [[--cosmic_ray | Installation:-Option-List#--cosmic_ray]].

<a name="OPT__FLAG_LOHNER_FORM"></a>
* #### `OPT__FLAG_LOHNER_FORM` &ensp; (1=FLASH-ghost1, 2=FLASH-ghost2, 3=form-invariant-ghost1, 4=form-invariant-ghost2) &ensp; [2]
    * **Description:**
Different forms of the Loehner's error estimator. This method
calculates the second derivative normalized by the first derivative, with
an additional filter term in the denominator to ignore small fluctuations.
`OPT__FLAG_LOHNER_FORM=2` corresponds to the form used by FLASH and
Enzo (see Eq. 8.4 in the [FLASH User Guide](http://flash.uchicago.edu/site/flashcode/user_support/flash4_ug_4p5/node59.html#SECTION05163000000000000000)
or Eq. 33 in the [Enzo Code Paper](http://iopscience.iop.org/article/10.1088/0067-0049/211/2/19/meta)).
`OPT__FLAG_LOHNER_FORM=1` is similar to `OPT__FLAG_LOHNER_FORM=2`
except using only 1 instead of 2 ghost zones on each side when
estimating spatial derivatives. `OPT__FLAG_LOHNER_FORM=3/4`
slightly revises the original formula to become form-invariant and
uses 1 and 2 ghost zones on each side when estimating spatial
derivatives, respectively.

        Specify the refinement thresholds on different levels in the input
file `Input__Flag_Lohner` with the
[[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_Lohner__FLASH2`
(must rename it as `Input__Flag_Lohner` to actually use it).

        ```
        # Level  Threshold_Refine  Threshold_Derefine    Filter    Soften    MinDensity
              0              0.80                0.80      0.01      0.00          0.00
              1              0.80                0.80      0.01      0.00          0.00
              2              0.80                0.80      0.01      0.00          0.00
        ```

        Note that this input table takes 6 columns, where the 1st column
(Level) is useless and the 5th column (Soften) is deprecated.

        * Threshold: dimensionless refinement/derefinement thresholds
(recommended values: 0.3 ~ 0.8).
        * Filter: dimensionless filter term in the denominator to
ignore small fluctuations (recommended value: 0.01 ~ 0.1).
        * Soften: dimensional floor value of the denominator. Specifically,
we calculate error = sqrt( N/max(D,Soften) ), where N and D are the
numerator and denominator in the Lohner's formula, respectively.
This parameter is deprecated.
        * MinDensity: minimum gas mass density threshold. All Loehner's
refinement criteria (i.e., all `OPT__FLAG_LOHNER_*` options) are disabled
for cells with gas density lower than this threshold.

      *Caution: currently all `OPT__FLAG_LOHNER_*` options share the same
input file `Input__Flag_Lohner`.*
    * **Restriction:**

<a name="OPT__FLAG_USER"></a>
* #### `OPT__FLAG_USER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: user-defined. Edit `src/Refine/Flag_User.cpp`
or a problem-specific function (for the latter, see
[[Add Problem Specific Functionalities | Adding-New-Simulations#vi-add-problem-specific-functionalities]]).
Specify the refinement thresholds on different levels in the input file
`Input__Flag_User` with the
[[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_User`.
See also [OPT__FLAG_USER_NUM](#OPT__FLAG_USER_NUM).
    * **Restriction:**

<a name="OPT__FLAG_USER_NUM"></a>
* #### `OPT__FLAG_USER_NUM` &ensp; (&#8805;1) &ensp; [1]
    * **Description:**
Number of threshold values on each level in the file `Input__Flag_User` for [OPT__FLAG_USER](#OPT__FLAG_USER).
    * **Restriction:**

<a name="OPT__FLAG_REGION"></a>
* #### `OPT__FLAG_REGION` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Edit `src/Refine/Flag_Region.cpp` to specify the regions that are
*allowed* to be refined. Note that this option does not trigger
any refinement. Instead, it simply forbids refinement outside the
specified regions.
    * **Restriction:**

<a name="OPT__FLAG_ANGULAR"></a>
* #### `OPT__FLAG_ANGULAR` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: angular resolution with respect to the specified center
([FLAG_ANGULAR_CEN_X](#FLAG_ANGULAR_CEN_X), [FLAG_ANGULAR_CEN_Y](#FLAG_ANGULAR_CEN_Y),
[FLAG_ANGULAR_CEN_Z](#FLAG_ANGULAR_CEN_Z)). Cells located at a distance greater than
`AngRes_Max_R` are not allowed to exceed the angular resolution `AngRes_Max` (in radians).
Cells are refined if their angular resolution is lower than `AngRes_Min` (in radians). Set
`AngRes_Max < 0.0` or `AngRes_Min < 0.0` to disable the respective criterion.
Specify the refinement thresholds `AngRes_Max, AngRes_Min, AngRes_Max_R`
on different levels in the input file `Input__Flag_AngularResolution`
with the [[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_AngularResolution`.
    * **Restriction:**
`AngRes_Max` has higher priority over `AngRes_Min`.
It is generally recommended to set `AngRes_Max < 0.5*AngRes_Min`.

<a name="FLAG_ANGULAR_CEN_X"></a>
* #### `FLAG_ANGULAR_CEN_X` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [box center]
    * **Description:**
x-coordinate of the center for calculating the angular resolution for [OPT__FLAG_ANGULAR](#OPT__FLAG_ANGULAR)
refinement criterion.
    * **Restriction:**

<a name="FLAG_ANGULAR_CEN_Y"></a>
* #### `FLAG_ANGULAR_CEN_Y` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [box center]
    * **Description:**
See [FLAG_ANGULAR_CEN_X](#FLAG_ANGULAR_CEN_X).
    * **Restriction:**

<a name="FLAG_ANGULAR_CEN_Z"></a>
* #### `FLAG_ANGULAR_CEN_Z` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [box center]
    * **Description:**
See [FLAG_ANGULAR_CEN_X](#FLAG_ANGULAR_CEN_X).
    * **Restriction:**

<a name="OPT__FLAG_RADIAL"></a>
* #### `OPT__FLAG_RADIAL` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: distance to the specified center ([FLAG_RADIAL_CEN_X](#FLAG_RADIAL_CEN_X),
[FLAG_RADIAL_CEN_Y](#FLAG_RADIAL_CEN_Y), [FLAG_RADIAL_CEN_Z](#FLAG_RADIAL_CEN_Z)).
Cells are refined if they are located at a distance smaller than `Refine_Rad`.
The value of `Refine_Rad` can be specified in the input file `Input__Flag_RadialResolution`
with the [[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_RadialResolution`.
    * **Restriction:**

<a name="FLAG_RADIAL_CEN_X"></a>
* #### `FLAG_RADIAL_CEN_X` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [box center]
    * **Description:**
x-coordinate of the center for calculating the radial resolution for [OPT__FLAG_RADIAL](#OPT__FLAG_RADIAL)
refinement criterion.
    * **Restriction:**

<a name="FLAG_RADIAL_CEN_Y"></a>
* #### `FLAG_RADIAL_CEN_Y` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [box center]
    * **Description:**
See [FLAG_RADIAL_CEN_X](#FLAG_RADIAL_CEN_X).
    * **Restriction:**

<a name="FLAG_RADIAL_CEN_Z"></a>
* #### `FLAG_RADIAL_CEN_Z` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [box center]
    * **Description:**
See [FLAG_RADIAL_CEN_X](#FLAG_RADIAL_CEN_X).
    * **Restriction:**

<a name="OPT__FLAG_NPAR_PATCH"></a>
* #### `OPT__FLAG_NPAR_PATCH` &ensp; (0=off, 1=itself, 2=itself+siblings) &ensp; [0]
    * **Description:**
Refinement criterion: number of particles in a patch.
For `OPT__FLAG_NPAR_PATCH=1`, only patches with more particles than a
given threshold are flagged for refinement. For `OPT__FLAG_NPAR_PATCH=2`,
not only patches with excessive numbers of particles but also their 26
siblings are flagged for refinement. Specify the refinement thresholds
on different levels in the input file `Input__Flag_NParPatch` with the
[[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_NParPatch`.
Note that the flag buffers (i.e.,
[FLAG_BUFFER_SIZE](#FLAG_BUFFER_SIZE),
[FLAG_BUFFER_SIZE_MAXM1_LV](#FLAG_BUFFER_SIZE_MAXM1_LV), and
[FLAG_BUFFER_SIZE_MAXM2_LV](#FLAG_BUFFER_SIZE_MAXM2_LV))
do not apply to this refinement criterion.
    * **Restriction:**
Currently always includes tracer particles.

<a name="OPT__FLAG_NPAR_CELL"></a>
* #### `OPT__FLAG_NPAR_CELL` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: number of particles in a cell.
Specify the refinement thresholds on different levels in the input file
`Input__Flag_NParCell` with the [[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_NParCell`.
    * **Restriction:**
Currently always excludes tracer particles.

<a name="OPT__FLAG_PAR_MASS_CELL"></a>
* #### `OPT__FLAG_PAR_MASS_CELL` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Refinement criterion: total particle mass in a cell.
Specify the refinement thresholds on different levels in the input file
`Input__Flag_ParMassCell` with the [[specific format | Runtime-Parameters:-Input__Flag_*]].
An example file can be found at `example/input/Input__Flag_ParMassCell`.
    * **Restriction:**

<a name="OPT__NO_FLAG_NEAR_BOUNDARY"></a>
* #### `OPT__NO_FLAG_NEAR_BOUNDARY` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Disallow refinement near the boundaries.
    * **Restriction:**

<a name="OPT__PATCH_COUNT"></a>
* #### `OPT__PATCH_COUNT` &ensp; (0=off, 1=every root-level step, 2=every substep) &ensp; [1]
    * **Description:**
Record the number of patches on each level in the log file
[[Record__PatchCount | Simulation Logs:-Record__Patchcount]].
    * **Restriction:**

<a name="OPT__PARTICLE_COUNT"></a>
* #### `OPT__PARTICLE_COUNT` &ensp; (0=off, 1=every root-level step, 2=every substep) &ensp; [1]
    * **Description:**
Record the number of particles on each level in the log file
[[Record__ParticleCount | Simulation Logs:-Record__Particlecount]].
    * **Restriction:**

<a name="OPT__REUSE_MEMORY"></a>
* #### `OPT__REUSE_MEMORY` &ensp; (0=off, 1=on, 2=aggressive) &ensp; [2]
    * **Description:**
Reuse allocated patch memory to reduce memory fragmentation.
For `OPT__REUSE_MEMORY=1`, the code will still deallocate patch memory
when redistributing all patches for load balancing
(see [[LB_INPUT__WLI_MAX | Runtime-Parameters:-MPI-and-OpenMP#LB_INPUT__WLI_MAX]]).
In comparison, for `OPT__REUSE_MEMORY=2`, the code will not deallocate
patch memory during the entire simulation. Note that this option will
not preallocate any patches unless [OPT__MEMORY_POOL](#OPT__MEMORY_POOL)
is enabled.
    * **Restriction:**

<a name="OPT__MEMORY_POOL"></a>
* #### `OPT__MEMORY_POOL` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Preallocate patches as a memory pool to reduce memory fragmentation.
One must specify the numbers of patches to be preallocated in the
input file [[Input__MemoryPool | Runtime-Parameters:-Input__MemoryPool]]
(check the link for details).
    * **Restriction:**
Only applicable when adopting [OPT__REUSE_MEMORY](#OPT__REUSE_MEMORY)=1/2.


## Remarks

### Potential outside the isolated boundaries
When adopting the isolated boundary conditions for gravity (i.e.,
[[OPT__BC_POT | Runtime-Parameters:-Gravity#OPT__BC_POT]]=2), the ghost zones of
gravitational potential outside the simulation domain are currently
filled out by extrapolation.


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
