#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_PointMass
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by Init_ExtPotAuxArray_PointMass(), where
//                      UserArray[0] = x coordinate of the external potential center
//                      UserArray[1] = y ...
//                      UserArray[2] = z ..
//                      UserArray[3] = gravitational_constant*point_source_mass
//                3. Currently it does not support the soften length
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_PointMass( const double x, const double y, const double z, const double Time, const double UserArray[] )
{

   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
   const real   GM     = (real)UserArray[3];
   const real   dx     = (real)(x - Cen[0]);
   const real   dy     = (real)(y - Cen[1]);
   const real   dz     = (real)(z - Cen[2]);
   const real   _r     = 1.0/SQRT( dx*dx + dy*dy + dz*dz );

   return -GM*_r;

} // FUNCTION : ExtPot_PointMass



// =================================
// get the CPU/GPU function pointers
// =================================

#ifdef __CUDACC__
__device__
#endif
static ExtPot_t ExtPot_Ptr = ExtPot_PointMass;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_PointMass
// Description :  Return the function pointers to the CPU/GPU external potential routines
//
// Note        :  1. To enable this routine, link to the function pointers "SetCPU/GPUExtPot_Ptr"
//                   in a test problem initializer as follows:
//
//                      void SetCPUExtPot_PointMass( ExtPot_t &CPUExtPot_Ptr );
//                      # ifdef GPU
//                      void SetGPUExtPot_PointMass( ExtPot_t &GPUExtPot_Ptr );
//                      # endif
//
//                      ...
//
//                      SetCPUExtPot_Ptr = SetCPUExtPot_PointMass;
//                      # ifdef GPU
//                      SetGPUExtPot_Ptr = SetGPUExtPot_PointMass;
//                      # endif
//
//                   --> Then it will be invoked by Init_ExtAccPot()
//                2. Must obtain the CPU and GPU function pointers by separate routines
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
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtPotAuxArray_PointMass
// Description :  Set the auxiliary array ExtPot_AuxArray[] used by ExtPot_PointMass()
//
// Note        :  1. To adopt this routine, link to the function pointer "Init_ExtPotAuxArray_Ptr"
//                   in a test problem initializer as follows:
//
//                      void Init_ExtPotAuxArray_PointMass( double AuxArray[] );
//
//                      ...
//
//                      Init_ExtPotAuxArray_Ptr = Init_ExtPotAuxArray_PointMass;
//
//                   --> Then it will be invoked by Init_ExtAccPot()
//                2. AuxArray[] has the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 10)
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void Init_ExtPotAuxArray_PointMass( double AuxArray[] )
{

// example parameters
   const double M  = 1.0;
   const double GM = NEWTON_G*M;

   AuxArray[0] = 0.5*amr->BoxSize[0];  // x coordinate of the external potential center
   AuxArray[1] = 0.5*amr->BoxSize[1];  // y ...
   AuxArray[2] = 0.5*amr->BoxSize[2];  // z ...
   AuxArray[3] = GM;                   // gravitational_constant*point_source_mass

} // FUNCTION : Init_ExtPotAuxArray_PointMass
#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
