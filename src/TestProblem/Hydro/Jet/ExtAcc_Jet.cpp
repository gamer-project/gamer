#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#endif

#ifdef GRAVITY




// =================================
// I. Set an auxiliary array
// =================================

#ifndef __CUDACC__

extern double Jet_BgVel[3];
extern double Jet_HSE_D;
extern double Jet_HSE_M200;
extern double Jet_HSE_R200;
extern double Jet_HSE_C200;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtAccAuxArray_Jet
// Description :  Set the auxiliary array ExtAcc_AuxArray[] used by ExtAcc_Jet()
//
// Note        :  1. Invoked by Init_ExtAcc_Jet()
//                2. AuxArray[] has the size of EXT_ACC_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray : Array to be filled up
//                Time     : Target physical time
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void SetExtAccAuxArray_Jet( double AuxArray[], const double Time )
{

   const double c = Jet_HSE_C200;

   AuxArray[0] = amr->BoxCenter[0];                                    // [0-2]: cluster center
   AuxArray[1] = amr->BoxCenter[1] - Jet_HSE_D;
   AuxArray[2] = amr->BoxCenter[2];
   AuxArray[3] = -NEWTON_G*Jet_HSE_M200/( log(1.0+c) - c/(1.0+c) );    // -G*M200/( log(1+c) - c/(1+c) )
   AuxArray[4] = Jet_HSE_R200 / c;                                     // scale radius
   AuxArray[5] = Jet_BgVel[1];

} // FUNCTION : SetExtAccAuxArray_Jet
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external acceleration
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_Jet
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by SetExtAccAuxArray_Jet()
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Target physical time
//                UserArray : User-provided auxiliary array
//
// Return      :  External acceleration Acc[] at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void ExtAcc_Jet( real Acc[], const double x, const double y, const double z, const double Time,
                        const double UserArray[] )
{

   const double Cen[3]  = { UserArray[0], UserArray[1], UserArray[2] };    // cluster center
   const real   coeff   = (real)UserArray[3];                              // -G*M200/( log(1+c) - c/(1+c) )
   const real   r0      = (real)UserArray[4];                              // scale radius
   const real   Vy      = (real)UserArray[5];                              // jet infall velocity along y
   const real   dx      = (real)0.0;                                       // do not consider dx for now
   const real   dy      = (real)(y - Cen[1] - Vy*Time);
   const real   dz      = (real)0.0;                                       // do not consider dz for now
   const real   r       = SQRT( dx*dx + dy*dy + dz*dz );
   const real   s       = r/r0;
   const real   one     = (real)1.0;
   const real   _r3     = one/CUBE(r);
   const real   acc_r   = coeff*( LOG(one+s) - s/(one+s) )*_r3;

   Acc[0] = acc_r*dx;
   Acc[1] = acc_r*dy;
   Acc[2] = acc_r*dz;

} // FUNCTION : ExtAcc_Jet



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtAcc_t ExtAcc_Ptr = ExtAcc_Jet;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_Jet
// Description :  Return the function pointers of the CPU/GPU external acceleration routines
//
// Note        :  1. Invoked by Init_ExtAcc_Jet()
//
// Parameter   :  CPU/GPUExtAcc_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtAcc_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtAcc_Jet( ExtAcc_t &GPUExtAcc_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtAcc_Ptr, ExtAcc_Ptr, sizeof(ExtAcc_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtAcc_Jet( ExtAcc_t &CPUExtAcc_Ptr )
{
   CPUExtAcc_Ptr = ExtAcc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtAccAuxArray_Jet( double [], const double );
void SetCPUExtAcc_Jet( ExtAcc_t & );
#ifdef GPU
void SetGPUExtAcc_Jet( ExtAcc_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtAcc_Jet
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
void Init_ExtAcc_Jet()
{

   SetExtAccAuxArray_Jet( ExtAcc_AuxArray, Time[0] );

   SetCPUExtAcc_Jet( CPUExtAcc_Ptr );
#  ifdef GPU
   SetGPUExtAcc_Jet( GPUExtAcc_Ptr );
#  endif

} // FUNCTION : Init_ExtAcc_Jet

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
