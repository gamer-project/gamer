## Compilation Options

Related options:
[[--gravity | Installation:-Option-List#--gravity]], &nbsp;
[[--pot_scheme | Installation:-Option-List#--pot_scheme]], &nbsp;
[[--store_pot_ghost | Installation:-Option-List#--store_pot_ghost]], &nbsp;
[[--unsplit_gravity | Installation:-Option-List#--unsplit_gravity]], &nbsp;
[[--comoving | Installation:-Option-List#--comoving]] &nbsp;


## Runtime Parameters
[[Runtime parameters: Gravity | Runtime-Parameters:-Gravity]]

Other related parameters:
[[DT__GRAVITY | Runtime-Parameters:-Timestep#DT__GRAVITY]], &nbsp;
[[POT_GPU_NPGROUP | Runtime-Parameters:-GPU#POT_GPU_NPGROUP]] &nbsp;


## Remarks

### External Acceleration/Potential

#### Using Analytical Function:

Follow the steps below to add external acceleration with an analytical
function. Get familiar with the general procedure of
[[Adding Problem-specific Functionalities | Adding-New-Simulations#vi-add-problem-specific-functionalities]]
first and be aware that adding external acceleration/potential takes a few extra steps
in order to utilize GPUs.

1. Enable [[OPT__EXT_ACC | Runtime-Parameters:-Gravity#OPT__EXT_ACC]].
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
Enable [[OPT__EXT_POT | Runtime-Parameters:-Gravity#OPT__EXT_POT]]`=1` and then follow the procedure
above and replace `Acc` by `Pot`. A point-mass example can be found
at `src/SelfGravity/CPU_Poisson/CPU_ExtPot_PointMass.cpp`.

**Caution: wave dark matter simulations
(i.e., [[--model | Installation:-Option-List#--model]]=ELBDM)
does not support external acceleration. Just use external potential.**

#### Using Table:

External potential can also be set via interpolating a table loaded from disk.
Set [[OPT__EXT_POT | Runtime-Parameters:-Gravity#OPT__EXT_POT]]`=2` and all `EXT_POT_TABLE_*` parameters.
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
* [[Main page of Physics Modules | Physics-Modules]]
* [[Main page of Runtime Parameters | Runtime-Parameters]]
