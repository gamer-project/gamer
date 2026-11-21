#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#endif

#ifdef GRAVITY






// =================================
// I. Set an auxiliary array
// =================================

#ifndef __CUDACC__

extern double g0;


//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtAccAuxArray_MTI
// Description :  Set the auxiliary array ExtAcc_AuxArray[] used by ExtAcc_MTI()
//
// Note        :  1. Invoked by Init_ExtAcc_MTI()
//                2. AuxArray[] has the size of EXT_ACC_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray : Array to be filled up
//                Time     : Target physical time
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void SetExtAccAuxArray_MTI( double AuxArray[], const double Time )
{

   AuxArray[0] = g0;

} // FUNCTION : SetExtAccAuxArray_MTI
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external acceleration
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_MTI
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by SetExtAccAuxArray_MTI()
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Target physical time
//                UserArray : User-provided auxiliary array
//
// Return      :  External acceleration Acc[] at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void ExtAcc_MTI( real Acc[], const double x, const double y, const double z, const double Time,
                        const double UserArray[] )
{

   Acc[0] = 0.0;
   Acc[1] = 0.0;
   Acc[2] = -(real)UserArray[0];

} // FUNCTION : ExtAcc_MTI



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtAcc_t ExtAcc_Ptr = ExtAcc_MTI;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_MTI
// Description :  Return the function pointers of the CPU/GPU external acceleration routines
//
// Note        :  1. Invoked by Init_ExtAcc_MTI()
//
// Parameter   :  CPU/GPUExtAcc_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtAcc_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtAcc_MTI( ExtAcc_t &GPUExtAcc_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtAcc_Ptr, ExtAcc_Ptr, sizeof(ExtAcc_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtAcc_MTI( ExtAcc_t &CPUExtAcc_Ptr )
{
   CPUExtAcc_Ptr = ExtAcc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtAccAuxArray_MTI( double [], const double );
void SetCPUExtAcc_MTI( ExtAcc_t & );
#ifdef GPU
void SetGPUExtAcc_MTI( ExtAcc_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtAcc_MTI
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
void Init_ExtAcc_MTI()
{

   SetExtAccAuxArray_MTI( ExtAcc_AuxArray, Time[0] );

   SetCPUExtAcc_MTI( CPUExtAcc_Ptr );
#  ifdef GPU
   SetGPUExtAcc_MTI( GPUExtAcc_Ptr );
#  endif

} // FUNCTION : Init_ExtAcc_MTI

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
