#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY



/********************************************************
1. Point-mass external acceleration
   --> It can be regarded as a template for implementing
       other external acceleration

2. This file is shared by both CPU and GPU

   GPU_Gravity/CUPOT_ExtAcc_PointMass.cu -> CPU_Gravity/CPU_ExtAcc_PointMass.cpp

3. Three steps are required to implement external acceleration

   I.   Set an auxiliary array
        --> SetExtAccAuxArray_PointMass()

   II.  Specify external acceleration
        --> ExtAcc_PointMass()

   III. Set initialization functions
        --> SetGPUExtAcc_PointMass()
            SetCPUExtAcc_PointMass()
            Init_ExtAcc_PointMass()

4. The external acceleration major routine, ExtAcc_PointMass(),
   must be thread-safe and not use any global variable

5. Reference: https://github.com/gamer-project/gamer/wiki/Gravity#external-accelerationpotential
********************************************************/



// soften length implementation
#  define SOFTEN_PLUMMER
//#  define SOFTEN_RUFFERT



// =================================
// I. Set an auxiliary array
// =================================

#ifndef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtAccAuxArray_PointMass
// Description :  Set the auxiliary array ExtAcc_AuxArray[] used by ExtAcc_PointMass()
//
// Note        :  1. Invoked by Init_ExtAcc_PointMass()
//                2. AuxArray[] has the size of EXT_ACC_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void SetExtAccAuxArray_PointMass( double AuxArray[] )
{

// example parameters
   const double M   = 1.0;
   const double GM  = NEWTON_G*M;
   const double Eps = 0.0;

   AuxArray[0] = 0.5*amr->BoxSize[0];  // x coordinate of the external acceleration center
   AuxArray[1] = 0.5*amr->BoxSize[1];  // y ...
   AuxArray[2] = 0.5*amr->BoxSize[2];  // z ...
   AuxArray[3] = GM;                   // gravitational_constant*point_source_mass
   AuxArray[4] = Eps;                  // soften_length (<=0.0 --> disable)

} // FUNCTION : SetExtAccAuxArray_PointMass
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external acceleration
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_PointMass
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by SetExtAccAuxArray_PointMass(), where
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
//                Time      : Target physical time
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
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtAcc_t ExtAcc_Ptr = ExtAcc_PointMass;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_PointMass
// Description :  Return the function pointers of the CPU/GPU external acceleration routines
//
// Note        :  1. Invoked by Init_ExtAcc_PointMass()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtAcc_PointMass( ExtAcc_t &CPUExtAcc_Ptr, ExtAcc_t &GPUExtAcc_Ptr )
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
}

#else // #ifdef __CUDACC__

void SetCPUExtAcc_PointMass( ExtAcc_t &CPUExtAcc_Ptr )
{
   CPUExtAcc_Ptr = ExtAcc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtAccAuxArray_PointMass( double [] );
void SetCPUExtAcc_PointMass( ExtAcc_t & );
#ifdef GPU
void SetGPUExtAcc_PointMass( ExtAcc_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtAcc_PointMass
// Description :  Initialize external acceleration
//
// Note        :  1. Set an auxiliary array by invoking SetExtAccAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external acceleration major routines by invoking SetCPU/GPUExtAcc_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtAcc_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtAcc_PointMass()
{

   SetExtAccAuxArray_PointMass( ExtAcc_AuxArray );

   SetCPUExtAcc_PointMass( CPUExtAcc_Ptr );
#  ifdef GPU
   SetGPUExtAcc_PointMass( GPUExtAcc_Ptr );
#  endif

} // FUNCTION : Init_ExtAcc_PointMass

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
