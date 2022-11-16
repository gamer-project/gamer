#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#endif

#ifdef GRAVITY



// soften length implementation
#  define SOFTEN_PLUMMER
//#  define SOFTEN_RUFFERT



// =================================
// I. Set an auxiliary array
// =================================

#ifndef __CUDACC__

extern double Bondi_MassBH;
extern double Bondi_Soften_R;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtAccAuxArray_Bondi
// Description :  Set the auxiliary array ExtAcc_AuxArray[] used by ExtAcc_Bondi()
//
// Note        :  1. Invoked by Init_ExtAcc_Bondi()
//                2. AuxArray[] has the size of EXT_ACC_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray : Array to be filled up
//                Time     : Target physical time
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void SetExtAccAuxArray_Bondi( double AuxArray[], const double Time )
{

   AuxArray[0] = amr->BoxCenter[0];
   AuxArray[1] = amr->BoxCenter[1];
   AuxArray[2] = amr->BoxCenter[2];
   AuxArray[3] = NEWTON_G*Bondi_MassBH;   // gravitational_constant*point_source_mass
   AuxArray[4] = Bondi_Soften_R;          // soften_length (<=0.0 --> disable)

} // FUNCTION : SetExtAccAuxArray_Bondi
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external acceleration
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_Bondi
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by SetExtAccAuxArray_Bondi()
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Target physical time
//                UserArray : User-provided auxiliary array
//
// Return      :  External acceleration Acc[] at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void ExtAcc_Bondi( real Acc[], const double x, const double y, const double z, const double Time,
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

} // FUNCTION : ExtAcc_Bondi



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtAcc_t ExtAcc_Ptr = ExtAcc_Bondi;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_Bondi
// Description :  Return the function pointers of the CPU/GPU external acceleration routines
//
// Note        :  1. Invoked by Init_ExtAcc_Bondi()
//
// Parameter   :  CPU/GPUExtAcc_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtAcc_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtAcc_Bondi( ExtAcc_t &GPUExtAcc_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtAcc_Ptr, ExtAcc_Ptr, sizeof(ExtAcc_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtAcc_Bondi( ExtAcc_t &CPUExtAcc_Ptr )
{
   CPUExtAcc_Ptr = ExtAcc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtAccAuxArray_Bondi( double [], const double );
void SetCPUExtAcc_Bondi( ExtAcc_t & );
#ifdef GPU
void SetGPUExtAcc_Bondi( ExtAcc_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtAcc_Bondi
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
void Init_ExtAcc_Bondi()
{

   SetExtAccAuxArray_Bondi( ExtAcc_AuxArray, Time[0] );

   SetCPUExtAcc_Bondi( CPUExtAcc_Ptr );
#  ifdef GPU
   SetGPUExtAcc_Bondi( GPUExtAcc_Ptr );
#  endif

} // FUNCTION : Init_ExtAcc_Bondi

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
