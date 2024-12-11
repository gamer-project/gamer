This page covers the following topics about constructing a simulation
initial condition (IC) and GAMER initialization:

* [Setting IC from Analytical Functions](#setting-ic-from-analytical-functions)
    * [Grids](#IC-Func-Grids)
    * [Magnetic Field](#IC-Func-BField)
    * [Particles](#IC-Func-Particles)
* [Setting IC from Files](#setting-ic-from-files)
    * [Grids](#IC-File-Grids)
    * [Magnetic Field](#IC-File-BField)
    * [Particles](#IC-File-Particles)
* [Compilation Options](#compilation-options)
* [Runtime Parameters](#runtime-parameters)


## Setting IC from Analytical Functions

<a name="IC-Func-Grids"></a>
### Grids

Set [OPT__INIT](#OPT__INIT)=1 and edit the following grid IC function:

* [[TESTPROB_ID | Runtime Parameters:-General#TESTPROB_ID]]=0:
edit the function `Init_Function_User()` in
`src/Model_Hydro/Hydro_Init_ByFunction_AssignData.cpp`.

* [[TESTPROB_ID | Runtime Parameters:-General#TESTPROB_ID]]&#8800;0:
edit a problem-specific initialization function
(usually named `SetGridIC()`). See also [[Adding New Simulations]].

The grid IC function (either `Init_Function_User()` or `SetGridIC()`)
has the following prototype:

```C++
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
```

It should set the variable `fluid[]` at a given location `x/y/z`
and time `Time`, where `fluid[]` is a 1D array to store different
fluid fields accessible by the keys `DENS`, `MOMX`, `MOMY`, `MOMZ`,
`ENGY` corresponding to gas mass density, momentum density along x/y/z,
and gas energy density (excluding potential energy), respectively.
The following example shows `SetGridIC()` of the blast wave test
(`src/TestProblem/Hydro/BlastWave/Init_TestProb_Hydro_BlastWave.cpp`):

```C++
//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid fields to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double Blast_Engy_Exp_Density = Blast_Engy_Exp/(4.0*M_PI/3.0*Blast_Radius*Blast_Radius*Blast_Radius);
   const double r = SQRT( SQR(x-Blast_Center[0]) + SQR(y-Blast_Center[1]) + SQR(z-Blast_Center[2]) );

   fluid[DENS] = Blast_Dens_Bg;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;

   if ( r <= Blast_Radius )   fluid[ENGY] = Blast_Engy_Exp_Density;
   else                       fluid[ENGY] = Blast_Engy_Bg;

} // FUNCTION : SetGridIC
```

If the IC does not have an analytical form and relies on the table
interpolation, one can use the helper functions
`Aux_LoadTable()` to load text tables from the disk and
`Mis_InterpolateFromTable()` to conduct linear interpolation.
See `src/TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp`
for example.

See also
[[ Add Problem-specific Grid Fields and Particle Attributes | Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes ]] for adding user-defined fluid fields.

> [!CAUTION]
> When enabling [[OPENMP | Installation:-Simulation Options#OPENMP]],
the grid IC function must be **thread-safe** since it will be invoked by
multiple threads in parallel.
>
> One can disable OpenMP parallelization
for the grid IC function by adopting
[[OPT__INIT_GRID_WITH_OMP| MPI-and-OpenMP#OPT__INIT_GRID_WITH_OMP]]=0.


<a name="IC-Func-BField"></a>
### Magnetic Field

One can either set the magnetic field directly or set the vector potential;
for the latter, the built-in routines will construct a corresponding
divergence-free magnetic field automatically. We describe the two approaches
separately in the following.

To set the magnetic field directly,
set [OPT__INIT](#OPT__INIT)=1 and [OPT__INIT_BFIELD_BYVECPOT](#OPT__INIT_BFIELD_BYVECPOT)=0.
Edit the following magnetic field IC function:

* [[TESTPROB_ID | Runtime Parameters:-General#TESTPROB_ID]]=0:
edit the function `Init_Function_BField_User()` in
`src/Model_Hydro/Hydro_Init_ByFunction_AssignData.cpp`.

* [[TESTPROB_ID | Runtime Parameters:-General#TESTPROB_ID]]&#8800;0:
edit a problem-specific initialization function
(usually named `SetBFieldIC()`). See also [[Adding New Simulations]].

The magnetic field IC function (either `Init_Function_BField_User()` or `SetBFieldIC()`)
has the following prototype:

```C++
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
```

It is similar to the grid IC function `SetGridIC()`. One should set the
variable `magnetic[]` at a given location `x/y/z` and time `Time`, where
`magnetic[]` is a 1D array to store different magnetic field components
accessible by the keys `MAGX`, `MAGY`, `MAGZ`.

To set the vector potential,
set [OPT__INIT](#OPT__INIT)=1 and [OPT__INIT_BFIELD_BYVECPOT](#OPT__INIT_BFIELD_BYVECPOT)=2.
Then edit a problem-specific initialization function linked to the function pointer
`Init_BField_ByVecPot_User_Ptr`, which has the following prototype:

```C++
double Init_BField_ByVecPot_User_Template( const double x, const double y, const double z, const double Time,
                                           const int lv, const char Component, double AuxArray[] )
```

It should return the vector potential at a given location `x/y/z` and time `Time`, where
`Component` specifies the vector potential component (`'x'`, `'y'`, or `'z'`) to be returned.

**Example: `Init_BField_ByVecPot_User_Template()` in `src/Model_Hydro/MHD_Init_BField_ByVecPot_Function.cpp`.**

> [!CAUTION]
> When enabling [[OPENMP | Installation:-Simulation Options#OPENMP]],
the magnetic field or the vector potential IC function must be **thread-safe** since it will be invoked by
multiple threads in parallel.
>
> One can disable OpenMP parallelization
for the magnetic field IC function by adopting
[[OPT__INIT_GRID_WITH_OMP| MPI-and-OpenMP#OPT__INIT_GRID_WITH_OMP]]=0.


<a name="IC-Func-Particles"></a>
### Particles

Set [[PAR_INIT | Particles#PAR_INIT]]=1 and edit the following
particle IC function:

* [[TESTPROB_ID | Runtime Parameters:-General#TESTPROB_ID]]=0:
edit the function `Par_Init_ByFunction()` in
`src/Particle/Par_Init_ByFunction.cpp`.

* [[TESTPROB_ID | Runtime Parameters:-General#TESTPROB_ID]]&#8800;0:
edit a problem-specific initialization function (e.g.,
`Par_Init_ByFunction_Merger()` in
`src/TestProblem/Hydro/ClusterMerger_vs_Flash/Par_Init_ByFunction_Merger.cpp`).
See also [[Adding New Simulations]].

The particle IC function has the following prototype:

```C++
void Par_Init_ByFunction( const long NPar_ThisRank, const long NPar_AllRank,
                          real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                          real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                          real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
```
It should set the particle IC in the arrays `ParMass`, `ParPosX/Y/Z`,
`ParVelX/Y/Z`, `ParTime`, `ParType`, and, optionally, the pointer array `*AllAttribute[PAR_NATT_TOTAL]`,
all of which have the size of `NPar_ThisRank` &#8212; the number of particles
to be set by this MPI rank. Note that particles set by this function
are only temporarily stored in this MPI rank and will later be
redistributed automatically to their host leaf patches. Therefore,
**there is no constraint on which particles should be set by this
function as long as each MPI rank initializes exactly `NPar_ThisRank`
particles**.

The built-in particle types (defined in `include/Macro.h`) include
`PTYPE_TRACER`, `PTYPE_GENERIC_MASSIVE`, `PTYPE_DARK_MATTER`, and `PTYPE_STAR`.
For `PTYPE_TRACER`, one must also enable the compilation option
[[TRACER | Installation: Simulation-Options#TRACER]].

The following example shows `Par_Init_ByFunction()` in
`src/Particle/Par_Init_ByFunction.cpp`:

```C++
//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  User-specified function to initialize particle attributes
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_LB_Init_RedistributeByRectangular()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                ParType       : Particle type     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data

//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction( const long NPar_ThisRank, const long NPar_AllRank,
                          real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                          real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                          real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{

// synchronize all particles to the physical time on the base level
// and assign particle type
   for (long p=0; p<NPar_ThisRank; p++)
   {
      ParTime[p] = Time[0];
      ParType[p] = PTYPE_GENERIC_MASSIVE;
   }

// set other particle attributes randomly
   real *ParPos[3] = { ParPosX, ParPosY, ParPosZ };
   real *ParVel[3] = { ParVelX, ParVelY, ParVelZ };

   const uint RSeed     = 2;                                         // random seed
   const real MassMin   = 1.0e-2;                                    // minimum value of particle mass
   const real MassMax   = 1.0;                                       // maximum value of particle mass
   const real PosMin[3] = { 0.0, 0.0, 0.0 };                         // minimum value of particle position
   const real PosMax[3] = { real( amr->BoxSize[0]*(1.0-1.0e-5) ),    // maximum value of particle position
                            real( amr->BoxSize[1]*(1.0-1.0e-5) ),
                            real( amr->BoxSize[2]*(1.0-1.0e-5) ) };
   const real VelMin[3] = { -1.0, -1.0, -1.0 };                      // minimum value of particle velocity
   const real VelMax[3] = { +1.0, +1.0, +1.0 };                      // maximum value of particle velocity

   srand( RSeed );

   for (long p=0; p<NPar_ThisRank; p++)
   {
      ParMass[p] = ( (real)rand()/RAND_MAX )*( MassMax - MassMin ) + MassMin;

      for (int d=0; d<3; d++)
      {
         ParPos[d][p] = ( (real)rand()/RAND_MAX )*( PosMax[d] - PosMin[d] ) + PosMin[d];
         ParVel[d][p] = ( (real)rand()/RAND_MAX )*( VelMax[d] - VelMin[d] ) + VelMin[d];
      }
   }

} // FUNCTION : Par_Init_ByFunction
```

See also
[[ Add Problem-specific Grid Fields and Particle Attributes | Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes ]] for adding user-defined particle attributes.

GAMER has a built-in routine for constructing a particle initial condition in equilibrium
(e.g., Plummer, NFW, tabular). Please refer to the test problem `Hydro/ParticleEquilibriumIC`
for its usage.


## Setting IC from Files

<a name="IC-File-Grids"></a>
### Grids
Set [OPT__INIT](#OPT__INIT)=3 to load the grid initial condition from a
uniform-mesh binary file named **`UM_IC`**. This file will be used to provide the
initial grid data of the entire computational domain _fully refined_
to the AMR level [OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL)
(but see also [OPT__UM_IC_DOWNGRADE](#OPT__UM_IC_DOWNGRADE) and
[OPT__UM_IC_REFINE](#OPT__UM_IC_REFINE) described below). The dimension of
this uniform-mesh file (assuming a row-major array) can be either
`[NFIELD][NZ][NY][NX]` (for [OPT__UM_IC_FORMAT](#OPT__UM_IC_FORMAT)=1) or
`[NZ][NY][NX][NFIELD]` (for [OPT__UM_IC_FORMAT](#OPT__UM_IC_FORMAT)=2), where
`NFIELD` is the number of fluid fields set by [OPT__UM_IC_NVAR](#OPT__UM_IC_NVAR)
and `NX/Y/Z` are the grid dimensions (i.e., number of cells) along the x/y/z
directions, respectively. Since `UM_IC` should store the initial condition of
a uniform mesh corresponding to level [OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL),
one must ensure that

* `NX`=[[NX0_TOT_X | Runtime-Parameters:-General#NX0_TOT_X]]*2^[OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL)
* `NY`=[[NX0_TOT_Y | Runtime-Parameters:-General#NX0_TOT_Y]]*2^[OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL)
* `NZ`=[[NX0_TOT_Z | Runtime-Parameters:-General#NX0_TOT_Z]]*2^[OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL)

For example, for `NX0_TOT_X=16`, `NX0_TOT_Y=32`, `NX0_TOT_Z=48`, `OPT__UM_IC_LEVEL=1`, and
`NCOMP_PASSIVE_USER=0`, `UM_IC` should have the dimension
`[5+0][48*2^1][32*2^1][16*2^1]=[5][96][64][32]`
for [OPT__UM_IC_FORMAT](#OPT__UM_IC_FORMAT)=1 or
`[96][64][32][5]` for [OPT__UM_IC_FORMAT](#OPT__UM_IC_FORMAT)=2,
assuming the array indices are row-major. The following C++ example
sets up a static and uniform gas with mass density of 1
and total energy density of 2 (assuming [OPT__UM_IC_FORMAT](#OPT__UM_IC_FORMAT)=1).

```c++
#include <cstdio>

int main()
{
   const int NX     = 32;
   const int NY     = 64;
   const int NZ     = 96;
   const int NFIELD = 5;

   float (*IC)[NZ][NY][NX] = new float [NFIELD][NZ][NY][NX];

   for (int k=0; k<NZ; k++)
   for (int j=0; j<NY; j++)
   for (int i=0; i<NX; i++)
   {
      IC[0][k][j][i] = 1.0;   // mass density
      IC[1][k][j][i] = 0.0;   // momentum density x
      IC[2][k][j][i] = 0.0;   // momentum density y
      IC[3][k][j][i] = 0.0;   // momentum density z
      IC[4][k][j][i] = 2.0;   // total energy density
   }

   FILE *File = fopen( "UM_IC", "wb" );
   fwrite( IC, sizeof(float), NFIELD*NZ*NY*NX, File );
   fclose( File );

   delete [] IC;
}
```

The entire computational domain will always be first fully refined to the
AMR level [OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL). After that, if
[OPT__UM_IC_DOWNGRADE](#OPT__UM_IC_DOWNGRADE) is enabled, the initialization
routine will remove grids on levels 1 &#8212; [OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL)
that do not satisfy any refinement criterion.
Also, if [OPT__UM_IC_REFINE](#OPT__UM_IC_REFINE) is enabled, the initialization
routine will add grids to levels [OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL)+1 &#8212;
[[ MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]
in the regions satisfying any refinement criterion.

The initial condition file `UM_IC` can be loaded concurrently by multiple
MPI processes using [OPT__UM_IC_LOAD_NRANK](#OPT__UM_IC_LOAD_NRANK).

To support AMR data in `UM_IC`, set [OPT__UM_IC_NLEVEL](#OPT__UM_IC_NLEVEL)>1 and
edit the input table `Input__UM_IC_RefineRegion`. See `example/input/Input__UM_IC_RefineRegion`
and the example code `tool/inits/create_UM_IC.cpp` for details.

> [!CAUTION]
> [OPT__INIT](#OPT__INIT)=3 does not fully support
user-defined passively advected scalars (i.e.,
[[NCOMP_PASSIVE_USER | Installation: Simulation-Options#NCOMP_PASSIVE_USER]]>0) yet.
[[Ask developers for help | Home#need-helps]] if needed.


<a name="IC-File-BField"></a>
### Magnetic Field
Set [OPT__INIT](#OPT__INIT)=1 and [OPT__INIT_BFIELD_BYVECPOT](#OPT__INIT_BFIELD_BYVECPOT)=1
to load the vector potential IC from a uniform-mesh HDF5 file
named **`B_IC`**. The built-in routines will construct a corresponding
divergence-free magnetic field automatically.

**Example: `tool/inits/gen_vec_pot.py`.**


<a name="IC-File-Particles"></a>
### Particles

Set [[PAR_INIT | Particles#PAR_INIT]]=3 to load the particle initial condition
from a binary file named **`PAR_IC`**. The dimension of this file (assuming a
row-major array) can be either
`[NUM_ATTRIBUTE][NUM_PARTICLE]` (for [[PAR_IC_FORMAT | Particles#par_ic_format]]=1) or
`[NUM_PARTICLE][NUM_ATTRIBUTE]` (for [[PAR_IC_FORMAT | Particles#par_ic_format]]=2), where
`NUM_ATTRIBUTE` is the number of particle attributes to be loaded
and `NUM_PARTICLE` is the total number of particles
(i.e., [[PAR_NPAR | Particles#par_npar]]).
By default, `NUM_ATTRIBUTE` is equal to
`8` + [[PAR_NATT_USER | Installation: Simulation-Options#PAR_NATT_USER]],
corresponding to particle mass, position x/y/z, velocity x/y/z,
type, and user-specified attributes (and in exactly this order).
One can also use [[PAR_IC_MASS | Particles#par_ic_mass]] / [[PAR_IC_TYPE | Particles#par_ic_type]]
to assign the same particle mass / type to all particles,
in which case the file `PAR_IC` should not store particle mass / type.

The following C++ example constructs a particle initial condition
file with 1000 particles assuming [[PAR_IC_MASS | Particles#par_ic_mass]]<0,
[[PAR_IC_TYPE | Particles#par_ic_type]]<0,
and [[PAR_IC_FORMAT | Particles#par_ic_format]]=1.

```c++
#include <cstdio>

int main()
{
   const int NUM_PARTICLE  = 1000;
   const int NUM_ATTRIBUTE = 8;

   float (*ParIC)[NUM_PARTICLE] = new float [NUM_ATTRIBUTE][NUM_PARTICLE];

   for (int p=0; p<NUM_PARTICLE; p++)
   {
//    replace the following lines by your particle initial condition
      ParIC[0][p] = 1.1;   // mass
      ParIC[1][p] = 2.2;   // position x
      ParIC[2][p] = 3.3;   // position y
      ParIC[3][p] = 4.4;   // position z
      ParIC[4][p] = 5.5;   // velocity x
      ParIC[5][p] = 6.6;   // velocity y
      ParIC[6][p] = 7.7;   // velocity z
      ParIC[7][p] = 1.0;   // type (generic massive)
   }

   FILE *File = fopen( "PAR_IC", "wb" );
   fwrite( ParIC, sizeof(float), NUM_PARTICLE*NUM_ATTRIBUTE, File );
   fclose( File );

   delete [] ParIC;
}
```

The built-in particle types (defined in `include/Macro.h`) include
`PTYPE_TRACER=0.0`, `PTYPE_GENERIC_MASSIVE=1.0`, `PTYPE_DARK_MATTER=2.0`, and `PTYPE_STAR=3.0`.
They have the floating-point types for now but will be changed to integers in the future.
For `PTYPE_TRACER`, one must also enable the compilation option
[[TRACER | Installation: Simulation-Options#TRACER]].

Note that it is not required to adopt [OPT__INIT](#OPT__INIT)=3 and
[[PAR_INIT | Particles#PAR_INIT]]=3 at the same time. In other words,
it is perfectly fine to set the grid initial condition from an analytical
function and load the particle initial condition from a file (and vice versa).


## Compilation Options

Related options:
[[NCOMP_PASSIVE_USER | Installation: Simulation-Options#NCOMP_PASSIVE_USER]], &nbsp;
[[PAR_NATT_USER | Installation: Simulation-Options#PAR_NATT_USER]] &nbsp;


## Runtime Parameters

Parameters described on this page:
[OPT__INIT](#OPT__INIT), &nbsp;
[OPT__INIT_BFIELD_BYVECPOT](#OPT__INIT_BFIELD_BYVECPOT), &nbsp;
[RESTART_LOAD_NRANK](#RESTART_LOAD_NRANK), &nbsp;
[OPT__RESTART_RESET](#OPT__RESTART_RESET), &nbsp;
[OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL), &nbsp;
[OPT__UM_IC_NLEVEL](#OPT__UM_IC_NLEVEL), &nbsp;
[OPT__UM_IC_NVAR](#OPT__UM_IC_NVAR), &nbsp;
[OPT__UM_IC_FORMAT](#OPT__UM_IC_FORMAT), &nbsp;
[OPT__UM_IC_FLOAT8](#OPT__UM_IC_FLOAT8), &nbsp;
[OPT__UM_IC_DOWNGRADE](#OPT__UM_IC_DOWNGRADE), &nbsp;
[OPT__UM_IC_REFINE](#OPT__UM_IC_REFINE), &nbsp;
[OPT__UM_IC_LOAD_NRANK](#OPT__UM_IC_LOAD_NRANK), &nbsp;
[OPT__INIT_RESTRICT](#OPT__INIT_RESTRICT), &nbsp;
[INIT_SUBSAMPLING_NCELL](#INIT_SUBSAMPLING_NCELL), &nbsp;
[OPT__FFTW_STARTUP](#OPT__FFTW_STARTUP) &nbsp;

Other related parameters:
[[PAR_INIT | Particles#PAR_INIT]], &nbsp;
[[PAR_IC_FORMAT | Particles#PAR_IC_FORMAT]], &nbsp;
[[PAR_IC_MASS | Particles#PAR_IC_MASS]], &nbsp;
[[OPT__INIT_GRID_WITH_OMP | MPI and OpenMP#OPT__INIT_GRID_WITH_OMP]] &nbsp;

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="OPT__INIT"></a>
* #### `OPT__INIT` &ensp; (1=function, 2=restart, 3=`UM_IC` file) &ensp; [none]
    * **Description:**
Grid initialization method.
`OPT__INIT=1`: using analytical functions; see
[Setting IC from Analytical Functions &#8212; Grids](#IC-Func-Grids).
`OPT__INIT=2`:
[[restarting from a simulation snapshot | Running-the-Code#restarting-from-a-snapshot]].
`OPT__INIT=3`: loading a uniform-mesh binary file named `UM_IC`; see
[Setting IC from Files &#8212; Grids](#IC-File-Grids).
    * **Restriction:**

<a name="OPT__INIT_BFIELD_BYVECPOT"></a>
* #### `OPT__INIT_BFIELD_BYVECPOT` &ensp; (0=off, 1=file, 2=function) &ensp; [0]
    * **Description:**
Set the magnetic field from either a vector potential file named `B_IC` (see
[Setting IC from Files &#8212; Magnetic Field](#IC-File-BField)) or an
analytical vector potential function (see
[Setting IC from Functions &#8212; Magnetic Field](#IC-Func-BField)).
    * **Restriction:**
For [[MHD | Installation: Simulation-Options#MHD]] only.

<a name="RESTART_LOAD_NRANK"></a>
* #### `RESTART_LOAD_NRANK` &ensp; (>0) &ensp; [1]
    * **Description:**
Number of parallel I/O for restart. In other words, `RESTART_LOAD_NRANK`
MPI processes will load the restart file in parallel.
    * **Restriction:**

<a name="OPT__RESTART_RESET"></a>
* #### `OPT__RESTART_RESET` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
During restart, reset some of the simulation status parameters
(e.g., step, time, snapshot ID) to their initial values as if the
simulation starts over again.
    * **Restriction:**

<a name="OPT__UM_IC_LEVEL"></a>
* #### `OPT__UM_IC_LEVEL` &ensp; (0 &#8804; input < [[NLEVEL | Installation:-Simulation Options#NLEVEL]]) &ensp; [0]
    * **Description:**
Starting AMR level in the uniform-mesh initial condition file.
See [Setting IC from Files &#8212; Grids](#IC-File-Grids) for details.
    * **Restriction:**

<a name="OPT__UM_IC_NLEVEL"></a>
* #### `OPT__UM_IC_NLEVEL` &ensp; (1 &#8804; input &#8804; [[NLEVEL | Installation:-Simulation Options#NLEVEL]]-[OPT__UM_IC_LEVEL](#OPT__UM_IC_LEVEL)) &ensp; [1]
    * **Description:**
Number of AMR levels in the uniform-mesh initial condition file.
See [Setting IC from Files &#8212; Grids](#IC-File-Grids) for details.
    * **Restriction:**

<a name="OPT__UM_IC_NVAR"></a>
* #### `OPT__UM_IC_NVAR` &ensp; (0 &#8804; input < total number of gas fields) &ensp; [total number of gas fields]
    * **Description:**
Number of fluid variables stored in the uniform-mesh initial condition file.
The default value is
5+[[NCOMP_PASSIVE_USER | Installation: Simulation-Options#NCOMP_PASSIVE_USER]]
for [[MODEL | Installation: Simulation-Options#MODEL]]=HYDRO.
See [Setting IC from Files &#8212; Grids](#IC-File-Grids) for details.
    * **Restriction:**

<a name="OPT__UM_IC_FORMAT"></a>
* #### `OPT__UM_IC_FORMAT` &ensp; (1=[v][z][y][x], 2=[z][y][x][v]; row-major and v=field) &ensp; [1]
    * **Description:**
Data format of the uniform-mesh initial condition file. See
[Setting IC from Files &#8212; Grids](#IC-File-Grids) for details.
    * **Restriction:**

<a name="OPT__UM_IC_FLOAT8"></a>
* #### `OPT__UM_IC_FLOAT8` &ensp; (<0: same as [[FLOAT8 | Installation:-Simulation Options#FLOAT8]], 0=single precision, 1=double precision) &ensp; [-1]
    * **Description:**
Floating-point precision of the uniform-mesh initial condition file.
    * **Restriction:**

<a name="OPT__UM_IC_DOWNGRADE"></a>
* #### `OPT__UM_IC_DOWNGRADE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Downgrade (i.e. derefine) the uniform-mesh initial condition data for cells
not satisfying any refinement criteria.
See [Setting IC from Files &#8212; Grids](#IC-File-Grids) for details.
    * **Restriction:**

<a name="OPT__UM_IC_REFINE"></a>
* #### `OPT__UM_IC_REFINE` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
Refine the uniform-mesh initial condition data from level `OPT__UM_IC_LEVEL` to
[[MAX_LEVEL | Runtime Parameters:-Refinement#MAX_LEVEL]] for cells satisfying the adopted
refinement criteria.
See [Setting IC from Files &#8212; Grids](#IC-File-Grids) for details.
    * **Restriction:**

<a name="OPT__UM_IC_LOAD_NRANK"></a>
* #### `OPT__UM_IC_LOAD_NRANK` &ensp; (>0) &ensp; [1]
    * **Description:**
Number of parallel I/O for loading the uniform-mesh initial condition file.
Specifically, it allows `OPT__UM_IC_LOAD_NRANK` MPI processes to load the
initial condition file concurrently. But the actually achieved parallel I/O
depends on the system specifications.
See also [Setting IC from Files &#8212; Grids](#IC-File-Grids).
    * **Restriction:**

<a name="OPT__INIT_RESTRICT"></a>
* #### `OPT__INIT_RESTRICT` &ensp; (0=off, 1=on) &ensp; [1]
    * **Description:**
For non-leaf patches, replace fluid data by the volume-weighted average
of their child patch data. It is similar to the option
[[OPT__FIXUP_RESTRICT | Hydro#OPT__FIXUP_RESTRICT]]
except that it only applies to the initial condition.
    * **Restriction:**

<a name="INIT_SUBSAMPLING_NCELL"></a>
* #### `INIT_SUBSAMPLING_NCELL` &ensp; (0=off, >0 &#8594; number of sub-cells along each direction) &ensp; [0]
    * **Description:**
Perform sub-sampling when constructing the grid IC to make it smoother.
Specifically, each cell will be divided into <var>N</var><sub>sub</sub><sup>3</sup>
sub-cells when calling the grid IC function, where
<var>N</var><sub>sub</sub> = `INIT_SUBSAMPLING_NCELL`, and then
take the volume-weighted average of these sub-cells.
    * **Restriction:**
Only applicable when adopting [OPT__INIT](#OPT__INIT)=1.

<a name="OPT__FFTW_STARTUP"></a>
* #### `OPT__FFTW_STARTUP` &ensp; (-1 &#8594; set to default, 0=ESTIMATE, 1=MEASURE, 2=PATIENT) &ensp; [-1]
    * **Description:**
Initialize FFTW plans. `MEASURE` is recommended for the balance
between FFTW plan initialization time and FFT performance.
Note that simulation results can vary in each run on the level of
machine precision for `OPT__FFTW_STARTUP != ESTIMATE`.
    * **Restriction:**
`PATIENT` is not supported by FFTW2.
Must use `ESTIMATE` when enabling
[[BITWISE_REPRODUCIBILITY | Installation: Simulation-Options#BITWISE_REPRODUCIBILITY]].


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Adding New Simulations | Adding New Simulations]]