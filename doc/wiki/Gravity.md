## Compilation Options

Related options:
[[--gravity | Installation:-Option-List#--gravity]], &nbsp;
[[--pot_scheme | Installation:-Option-List#--pot_scheme]], &nbsp;
[[--store_pot_ghost | Installation:-Option-List#--store_pot_ghost]], &nbsp;
[[--unsplit_gravity | Installation:-Option-List#--unsplit_gravity]], &nbsp;
[[--comoving | Installation:-Option-List#--comoving]] &nbsp;


## Runtime Parameters

Parameters described on this page:
[OPT__BC_POT](#OPT__BC_POT), &nbsp;
[GFUNC_COEFF0](#GFUNC_COEFF0), &nbsp;
[NEWTON_G](#NEWTON_G), &nbsp;
[SOR_OMEGA](#SOR_OMEGA), &nbsp;
[SOR_MAX_ITER](#SOR_MAX_ITER), &nbsp;
[SOR_MIN_ITER](#SOR_MIN_ITER), &nbsp;
[MG_MAX_ITER](#MG_MAX_ITER), &nbsp;
[MG_NPRE_SMOOTH](#MG_NPRE_SMOOTH), &nbsp;
[MG_NPOST_SMOOTH](#MG_NPOST_SMOOTH), &nbsp;
[MG_TOLERATED_ERROR](#MG_TOLERATED_ERROR), &nbsp;
[OPT__GRA_P5_GRADIENT](#OPT__GRA_P5_GRADIENT), &nbsp;
[OPT__SELF_GRAVITY](#OPT__SELF_GRAVITY), &nbsp;
[OPT__EXT_ACC](#OPT__EXT_ACC), &nbsp;
[OPT__EXT_POT](#OPT__EXT_POT), &nbsp;
[EXT_POT_TABLE_NAME](#EXT_POT_TABLE_NAME), &nbsp;
[EXT_POT_TABLE_NPOINT_X](#EXT_POT_TABLE_NPOINT_X), &nbsp;
[EXT_POT_TABLE_NPOINT_Y](#EXT_POT_TABLE_NPOINT_Y), &nbsp;
[EXT_POT_TABLE_NPOINT_Z](#EXT_POT_TABLE_NPOINT_Z), &nbsp;
[EXT_POT_TABLE_DH](#EXT_POT_TABLE_DH), &nbsp;
[EXT_POT_TABLE_EDGEL_X](#EXT_POT_TABLE_EDGEL_X), &nbsp;
[EXT_POT_TABLE_EDGEL_Y](#EXT_POT_TABLE_EDGEL_Y), &nbsp;
[EXT_POT_TABLE_EDGEL_Z](#EXT_POT_TABLE_EDGEL_Z)

Other related parameters:
[[DT__GRAVITY | Runtime-Parameters:-Timestep#DT__GRAVITY]], &nbsp;
[[POT_GPU_NPGROUP | GPU#POT_GPU_NPGROUP]] &nbsp;

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

 <a name="OPT__BC_POT"></a>
* #### `OPT__BC_POT` &ensp; (1=periodic, 2=isolated) &ensp; [none]
    * **Description:**
Gravity boundary condition. See also
[[Potential Outside the Isolated Boundaries | Runtime-Parameters:-Refinement#potential-outside-the-isolated-boundaries]].
    * **Restriction:**

<a name="GFUNC_COEFF0"></a>
* #### `GFUNC_COEFF0` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [depend]
    * **Description:**
Green's function coefficient for calculating the "self"
gravitational potential of each cell (i.e., the gravitational
potential contribution from the mass within the same cell).
The default value depends on the particle interpolation scheme
(`PAR_INTERP`). Set GFUNC_COEFF0=0.0 to exclude the self-potential.
    * **Restriction:**
Only applicable to the isolated boundary condition.

<a name="NEWTON_G"></a>
* #### `NEWTON_G` &ensp; (>0.0) &ensp; [conform to the unit system set by [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]]]
    * **Description:**
Gravitational constant.
    * **Restriction:**
It will be overwritten by the default value when [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]]
is on; no default when [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]] is off.

<a name="SOR_OMEGA"></a>
* #### `SOR_OMEGA` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [1.69]
    * **Description:**
Parameter of the SOR Poisson solver for optimizing the convergence rate.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=SOR.

<a name="SOR_MAX_ITER"></a>
* #### `SOR_MAX_ITER` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [single precision=60, double precision=110]
    * **Description:**
Maximum number of iterations in the SOR Poisson solver.
The default value depends on the adopted floating-point accuracy (`FLOAT8`).
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=SOR.

<a name="SOR_MIN_ITER"></a>
* #### `SOR_MIN_ITER` &ensp; (&#8805;3; <0 &#8594; set to default) &ensp; [10]
    * **Description:**
Minimum number of iterations in the SOR Poisson solver.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=SOR.

<a name="MG_MAX_ITER"></a>
* #### `MG_MAX_ITER` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [single precision=10, double precision=20]
    * **Description:**
Maximum number of iterations in the multigrid Poisson solver.
The default value depends on the adopted floating-point accuracy ([[--double | Installation:-Option-List#--double]]).
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=MG.

<a name="MG_NPRE_SMOOTH"></a>
* #### `MG_NPRE_SMOOTH` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [3]
    * **Description:**
Number of pre-smoothing steps in the multigrid Poisson solver.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=MG.

<a name="MG_NPOST_SMOOTH"></a>
* #### `MG_NPOST_SMOOTH` &ensp; (&#8805;0; <0 &#8594; set to default) &ensp; [3]
    * **Description:**
Number of post-smoothing steps in the multigrid Poisson solver.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=MG.

<a name="MG_TOLERATED_ERROR"></a>
* #### `MG_TOLERATED_ERROR` &ensp; (&#8805;0.0; <0.0 &#8594; set to default) &ensp; [single precision=1e-6, double precision=1e-15]
    * **Description:**
Maximum tolerable error in the multigrid Poisson solver.
    * **Restriction:**
Only applicable when adopting the compilation option
[[--pot_scheme | Installation:-Option-List#--pot_scheme]]=MG.

<a name="OPT__GRA_P5_GRADIENT"></a>
* #### `OPT__GRA_P5_GRADIENT` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Use 5-point instead of 3-point gradient for calculating the
gravitational acceleration from potential.
*This functionality is only experimental.*
    * **Restriction:**
Must manually set `#define GRA_GHOST_SIZE 2` (and `#define USG_GHOST_SIZE 2`
as well when adopting the compilation option
[[--unsplit_gravity | Installation:-Option-List#--unsplit_gravity]])
in the header file `Macro.h`. Unsupported for particle update.

<a name="OPT__SELF_GRAVITY"></a>
* #### `OPT__SELF_GRAVITY` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Enable self-gravity.
    * **Restriction:**

<a name="OPT__EXT_ACC"></a>
* #### `OPT__EXT_ACC` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Add external acceleration. See
[External Acceleration/Potential](#external-accelerationpotential)
for how to specify external acceleration.
    * **Restriction:**
Not applicable to the wave dark matter simulations
([[--model | Installation:-Option-List#--model]]=ELBDM).

<a name="OPT__EXT_POT"></a>
* #### `OPT__EXT_POT` &ensp; (0=off, 1=function, 2=table) &ensp; [0]
    * **Description:**
Add external potential. See
[External Acceleration/Potential](#external-accelerationpotential)
for how to specify external potential.
    * **Restriction:**

<a name="EXT_POT_TABLE_NAME"></a>
* #### `EXT_POT_TABLE_NAME` &ensp; (string) &ensp; [none]
    * **Description:**
For [OPT__EXT_POT](#OPT__EXT_POT)`=2`: table name.
    * **Restriction:**

<a name="EXT_POT_TABLE_NPOINT_X"></a>
* #### `EXT_POT_TABLE_NPOINT_X` &ensp; (>1) &ensp; [none]
    * **Description:**
For [OPT__EXT_POT](#OPT__EXT_POT)`=2`: table size (i.e., number of data points) along the x direction.
    * **Restriction:**

<a name="EXT_POT_TABLE_NPOINT_Y"></a>
* #### `EXT_POT_TABLE_NPOINT_Y` &ensp; (>1) &ensp; [none]
    * **Description:**
See [EXT_POT_TABLE_NPOINT_X](#EXT_POT_TABLE_NPOINT_X).
    * **Restriction:**

<a name="EXT_POT_TABLE_NPOINT_Z"></a>
* #### `EXT_POT_TABLE_NPOINT_Z` &ensp; (>1) &ensp; [none]
    * **Description:**
See [EXT_POT_TABLE_NPOINT_X](#EXT_POT_TABLE_NPOINT_X).
    * **Restriction:**

<a name="EXT_POT_TABLE_DH"></a>
* #### `EXT_POT_TABLE_DH` &ensp; (>0.0) &ensp; [none]
    * **Description:**
For [OPT__EXT_POT](#OPT__EXT_POT)`=2`: spatial interval between adjacent tabular data.
    * **Restriction:**

<a name="EXT_POT_TABLE_EDGEL_X"></a>
* #### `EXT_POT_TABLE_EDGEL_X` &ensp; (cover the simulation domain plus ghost zones) &ensp; [none]
    * **Description:**
For [OPT__EXT_POT](#OPT__EXT_POT)`=2`: starting x coordinate of the tabular data.
    * **Restriction:**
Must cover the simulation domain plus gravity ghost zones. Specifically, `EXT_POT_TABLE_EDGEL_X` must be
smaller or equal to `DomainLeftEdge - 1*RootLevelCell` and
`EXT_POT_TABLE_EDGEL_X + (EXT_POT_TABLE_NPOINT_X-1)*EXT_POT_TABLE_DH`
must be larger or equal to `DomainRightEdge + 1*RootLevelCell`.

<a name="EXT_POT_TABLE_EDGEL_Y"></a>
* #### `EXT_POT_TABLE_EDGEL_Y` &ensp; (cover the simulation domain) &ensp; [none]
    * **Description:**
See [EXT_POT_TABLE_EDGEL_X](#EXT_POT_TABLE_EDGEL_X).
    * **Restriction:**

<a name="EXT_POT_TABLE_EDGEL_Z"></a>
* #### `EXT_POT_TABLE_EDGEL_Z` &ensp; (cover the simulation domain) &ensp; [none]
    * **Description:**
See [EXT_POT_TABLE_EDGEL_X](#EXT_POT_TABLE_EDGEL_X).
    * **Restriction:**



## Remarks

### External Acceleration/Potential

#### Using Analytical Function:

Follow the steps below to add external acceleration with an analytical
function. Get familiar with the general procedure of
[[Adding Problem-specific Functionalities | Adding-New-Simulations#vi-add-problem-specific-functionalities]]
first and be aware that adding external acceleration/potential takes a few extra steps
in order to utilize GPUs.

1. Enable [OPT__EXT_ACC](#OPT__EXT_ACC).
2. Go to your new test problem folder and copy the built-in point-mass acceleration
routine as a template:
   ```
   cp ../../../SelfGravity/CPU_Gravity/CPU_ExtAcc_PointMass.cpp ExtAcc_NewProblem.cpp
   ```
3. Edit `ExtAcc_NewProblem.cpp` to
   * Replace the keyword `PointMass` by `NewProblem` (or whatever appropriate).
   * Edit the routine `SetExtAccAuxArray_NewProblem()` to set the
     auxiliary array `AuxArray[]`, which will be automatically
     passed to `ExtAcc_NewProblem()` as the input array `UserArray[]`.
     The array size is set by `EXT_ACC_NAUX_MAX` in `include/Macro.h` (default 20).
   * Edit the routine `ExtAcc_NewProblem()`, which should return the external
     acceleration `Acc[]` at `(x, y, z, Time)`.
4. Link `ExtAcc_NewProblem.cpp` to a corresponding CUDA source file:
   ```
   ln -s ExtAcc_NewProblem.cpp ExtAcc_NewProblem.cu
   ```
   This avoids writing redundant codes for CPU and GPU runs separately.
5. Edit your problem source file (e.g., `Init_TestProb_Hydro_NewProblem.cpp`) to
   * Declare a function prototype on top of the file:
      ```C++
      void Init_ExtAcc_NewProblem();
      ```
   * Set up a function pointer in your problem initialization
   function (e.g., `Init_TestProb_Hydro_NewProblem()`):
      ```C++
      # ifdef GRAVITY
        Init_ExtAcc_Ptr = Init_ExtAcc_NewProblem
      # endif
      ```

Adding external potential can be done in a very similar way.
Enable [OPT__EXT_POT](#OPT__EXT_POT)`=1` and then follow the procedure
above and replace `Acc` by `Pot`. A point-mass example can be found
at `src/SelfGravity/CPU_Poisson/CPU_ExtPot_PointMass.cpp`.

**Caution: wave dark matter simulations
(i.e., [[--model | Installation:-Option-List#--model]]=ELBDM)
does not support external acceleration. Just use external potential.**

#### Using Table:

External potential can also be set via interpolating a table loaded from disk.
Set [OPT__EXT_POT](#OPT__EXT_POT)`=2` and all `EXT_POT_TABLE_*` parameters.
The following C++ example creates such a table with a point-mass potential.

```C++
#include <cstdio>
#include <cmath>

typedef double real; // this must be consistent with the compilation option FLOAT8

int main()
{
// parameters
   const char  *Filename = "ExtPotTable";       // EXT_POT_TABLE_NAME
   const int    NX       = 131;                 // EXT_POT_TABLE_NPOINT_X
   const int    NY       = 131;                 // EXT_POT_TABLE_NPOINT_Y
   const int    NZ       = 131;                 // EXT_POT_TABLE_NPOINT_Z
   const double dh       = 0.0078125;           // EXT_POT_TABLE_NPOINT_DH
   const double EdgeL[3] = { -dh, -dh, -dh };   // EXT_POT_TABLE_EDGEL_X/Y/Z

// set the external potential array
   const double G  = 1.0;                       // gravitational constant
   const double M  = 1.0;                       // particle mass
   const double cx = EdgeL[0] + 0.5*(NX-1)*dh;  // centre
   const double cy = EdgeL[1] + 0.5*(NY-1)*dh;
   const double cz = EdgeL[2] + 0.5*(NZ-1)*dh;

   double dx, dy, dz, r;

   real (*ExtPot)[NY][NX] = new real [NZ][NY][NX];

   for (int k=0; k<NZ; k++)  {  dz = EdgeL[2] + k*dh - cz;
   for (int j=0; j<NY; j++)  {  dy = EdgeL[1] + j*dh - cy;
   for (int i=0; i<NX; i++)  {  dx = EdgeL[0] + i*dh - cx;

      r               = sqrt( dx*dx + dy*dy + dz*dz );
      ExtPot[k][j][i] = -G*M/r;

   }}}

// output
   FILE *File = fopen( Filename, "wb" );
   fwrite( ExtPot, sizeof(real), (long)NX*NY*NZ, File );
   fclose( File );

   delete [] ExtPot;
}
```

[Optional] To implement a user-specified interpolation routine, use
`src/SelfGravity/CPU_Poisson/CPU_ExtPot_Tabular.cpp` as a template and
follow the steps [2-5] in [Using Analytical Function](#using-analytical-function).

**Caution: external acceleration does not support this functionality yet.**


<br>

## Links
* [[Main page of Runtime Parameters | Runtime-Parameters]]
