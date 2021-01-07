#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY



/********************************************************
1. Point-mass external potential
   --> It can be regarded as a template for implementing
       other external potential

2. This file is shared by both CPU and GPU

   GPU_Poisson/CUPOT_ExtPot_PointMass.cu -> CPU_Poisson/CPU_ExtPot_PointMass.cpp

3. Three steps are required to implement external potential

   I.   Set auxiliary arrays
        --> SetExtPotAuxArray_PointMass()

   II.  Specify external potential
        --> ExtPot_PointMass()

   III. Set initialization functions
        --> SetGPUExtPot_PointMass()
            SetCPUExtPot_PointMass()
            Init_ExtPot_PointMass()

4. The external potential major routine, ExtPot_PointMass(),
   must be thread-safe and not use any global variable

5. Reference: https://github.com/gamer-project/gamer/wiki/Gravity#external-accelerationpotential
********************************************************/



// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_PointMass
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_PointMass()
//
// Note        :  1. Invoked by Init_ExtPot_PointMass()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_PointMass( double AuxArray_Flt[], int AuxArray_Int[] )
{

// example parameters
   const double M  = 1.0;
   const double GM = NEWTON_G*M;

   AuxArray_Flt[0] = 0.5*amr->BoxSize[0];    // x coordinate of the external potential center
   AuxArray_Flt[1] = 0.5*amr->BoxSize[1];    // y ...
   AuxArray_Flt[2] = 0.5*amr->BoxSize[2];    // z ...
   AuxArray_Flt[3] = GM;                     // gravitational_constant*point_source_mass

} // FUNCTION : SetExtPotAuxArray_PointMass
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_PointMass
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_PointMass(), where
//                      UserArray_Flt[0] = x coordinate of the external potential center
//                      UserArray_Flt[1] = y ...
//                      UserArray_Flt[2] = z ..
//                      UserArray_Flt[3] = gravitational_constant*point_source_mass
//                3. Currently it does not support the soften length
//
// Parameter   :  x/y/z             : Target spatial coordinates
//                Time              : Target physical time
//                UserArray_Flt/Int : User-provided floating-point/integer auxiliary arrays
//                Usage             : Different usages of external potential when computing total potential on level Lv
//                                    --> EXT_POT_USAGE_ADD     : add external potential on Lv
//                                        EXT_POT_USAGE_SUB     : subtract external potential for preparing self-gravity potential on Lv-1
//                                        EXT_POT_USAGE_SUB_TINT: like SUB but for temporal interpolation
//                                    --> This parameter is useless in most cases
//                PotTable          : 3D potential table used by EXT_POT_TABLE
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_PointMass( const double x, const double y, const double z, const double Time,
                              const double UserArray_Flt[], const int UserArray_Int[],
                              const ExtPotUsage_t Usage, const real PotTable[] )
{

   const double Cen[3] = { UserArray_Flt[0], UserArray_Flt[1], UserArray_Flt[2] };
   const real   GM     = (real)UserArray_Flt[3];
   const real   dx     = (real)(x - Cen[0]);
   const real   dy     = (real)(y - Cen[1]);
   const real   dz     = (real)(z - Cen[2]);
   const real   _r     = 1.0/SQRT( dx*dx + dy*dy + dz*dz );

   return -GM*_r;

} // FUNCTION : ExtPot_PointMass



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_PointMass;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_PointMass
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_PointMass()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtPot_PointMass( ExtPot_t &CPUExtPot_Ptr, ExtPot_t &GPUExtPot_Ptr )
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_PointMass( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_PointMass( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_PointMass( double [], int [] );
void SetCPUExtPot_PointMass( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_PointMass( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_PointMass
// Description :  Initialize external potential
//
// Note        :  1. Set auxiliary arrays by invoking SetExtPotAuxArray_*()
//                   --> They will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external potential major routines by invoking SetCPU/GPUExtPot_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtPot_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtPot_PointMass()
{

   SetExtPotAuxArray_PointMass( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int );
   SetCPUExtPot_PointMass( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_PointMass( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_PointMass

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
