#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY



// soften length implementation
#  define SOFTEN_PLUMMER
//#  define SOFTEN_RUFFERT




//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_PointMass
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by Init_ExtAccAuxArray_PointMass(), where
//                      UserArray[0] = x coordinate of the external acceleration center
//                      UserArray[1] = y ...
//                      UserArray[2] = z ..
//                      UserArray[3] = gravitational_constant*point_source_mass
//                      UserArray[4] = soften_length (<=0.0 --> disable)
//                3. Two different soften length implementations are supported
//                   --> SOFTEN_PLUMMER & SOFTEN_RUFFERT
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Current physical time
//                UserArray : User-provided auxiliary array
//
// Return      :  External acceleration Acc[] at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void ExtAcc_PointMass( real Acc[], const double x, const double y, const double z, const double Time,
                              const double UserArray[] )
{

   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
   const real GM       = (real)UserArray[3];
   const real eps      = (real)UserArray[4];
   const real dx       = (real)(x - Cen[0]);
   const real dy       = (real)(y - Cen[1]);
   const real dz       = (real)(z - Cen[2]);
   const real r        = SQRT( dx*dx + dy*dy + dz*dz );

// Plummer
#  if   ( defined SOFTEN_PLUMMER )
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps), (real)-1.5 );

// Ruffert 1994
#  elif ( defined SOFTEN_RUFFERT )
   const real tmp = EXP( -SQR(r)/SQR(eps) );
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps)*tmp, (real)-1.5 )*( (real)1.0 - tmp );

#  else
   const real _r3 = (real)1.0/CUBE(r);
#  endif

   Acc[0] = -GM*_r3*dx;
   Acc[1] = -GM*_r3*dy;
   Acc[2] = -GM*_r3*dz;

} // FUNCTION : ExtAcc_PointMass



// =================================
// get the CPU/GPU function pointers
// =================================

#ifdef __CUDACC__
__device__
#endif
static ExtAcc_t ExtAcc_Ptr = ExtAcc_PointMass;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_PointMass
// Description :  Return the function pointers to the CPU/GPU external acceleration routines
//
// Note        :  1. To enable this routine, link to the function pointers "SetCPU/GPUExtAcc_Ptr"
//                   in a test problem initializer as follows:
//
//                      void SetCPUExtAcc_PointMass( ExtAcc_t &CPUExtAcc_Ptr );
//                      # ifdef GPU
//                      void SetGPUExtAcc_PointMass( ExtAcc_t &GPUExtAcc_Ptr );
//                      # endif
//
//                      ...
//
//                      SetCPUExtAcc_Ptr = SetCPUExtAcc_PointMass;
//                      # ifdef GPU
//                      SetGPUExtAcc_Ptr = SetGPUExtAcc_PointMass;
//                      # endif
//
//                   --> Then it will be invoked by Init_ExtAccPot()
//
// Parameter   :  CPU/GPUExtAcc_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtAcc_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtAcc_PointMass( ExtAcc_t &GPUExtAcc_Ptr )
{

   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtAcc_Ptr, ExtAcc_Ptr, sizeof(ExtAcc_t) )  );

} // FUNCTION : GetGPUFuncPtr_ExtAcc

#else // #ifdef __CUDACC__

void SetCPUExtAcc_PointMass( ExtAcc_t &CPUExtAcc_Ptr )
{

   CPUExtAcc_Ptr = ExtAcc_Ptr;

} // FUNCTION : GetCPUFuncPtr_ExtAcc

#endif // #ifdef __CUDACC__ ... else ...



#endif // #ifdef GRAVITY
