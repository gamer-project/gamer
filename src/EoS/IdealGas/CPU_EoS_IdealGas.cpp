#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#if ( MODEL == HYDRO )



/**********************************************
This file is shared by both CPU and GPU

    CPU_EoS_IdealGas.cpp
    GPU_EoS_IdealGas.cu

Three steps are required to implement an EoS

   I.   Specify EoS
   II.  Set the CPU/GPU function pointers
   III. Set the auxiliary array
**********************************************/



// =====================================
// I. Specify EoS
// =====================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_IdealGas
// Description :  Convert gas mass density and internal energy density to gas pressure using an ideal-gas EoS
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. Auxiliary array UserArray[] is set by EoS_InitAuxArray_IdealGas(), where
//                      UserArray[0] = gamma
//                      UserArray[1] = gamma-1
//                      UserArray[2] = 1/(gamma-1)
//
// Parameter   :  Dens      : Gas mass density
//                Eint      : Gas internal energy density
//                UserArray : User-provided auxiliary array (see the Note above)
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Pres_IdealGas( const real Dens, const real Eint, const double UserArray[] )
{

   const real Gamma_m1 = (real)UserArray[1];
   real Pres;

   Pres = Eint * Gamma_m1;

   return Pres;

} // FUNCTION : EoS_DensEint2Pres_IdealGas



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2Eint_IdealGas
// Description :  Convert gas mass density and pressure to gas internal energy density using an ideal-gas EoS
//
// Note        :  1. See EoS_DensEint2Pres_IdealGas()
//
// Parameter   :  Dens      : Gas mass density
//                Pres      : Gas pressure
//                UserArray : User-provided auxiliary array (see the Note above)
//
// Return      :  Gas internal energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2Eint_IdealGas( const real Dens, const real Pres, const double UserArray[] )
{

   const real _Gamma_m1 = (real)UserArray[2];
   real Eint;

   Eint = Pres * _Gamma_m1;

   return Eint;

} // FUNCTION : EoS_DensPres2Eint_IdealGas



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_IdealGas
// Description :  Convert gas mass density and pressure to sound speed squared using an ideal-gas EoS
//
// Note        :  1. See EoS_DensEint2Pres_IdealGas()
//
// Parameter   :  Dens      : Gas mass density
//                Pres      : Gas pressure
//                UserArray : User-provided auxiliary array (see the Note above)
//
// Return      :  Sound speed square
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_IdealGas( const real Dens, const real Pres, const double UserArray[] )
{

   const real Gamma = (real)UserArray[0];
   real Cs2;

   Cs2 = Gamma * Pres / Dens;

   return Cs2;

} // FUNCTION : EoS_DensPres2CSqr_IdealGas



// =====================================
// II. Set the CPU/GPU function pointers
// =====================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_IdealGas;
FUNC_SPACE EoS_DP2E_t EoS_DensPres2Eint_Ptr = EoS_DensPres2Eint_IdealGas;
FUNC_SPACE EoS_DP2C_t EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_IdealGas;



//-----------------------------------------------------------------------------------------
// Function    :  EoS_InitCPU/GPUFunc_IdealGas
// Description :  Return the function pointers to the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init() when adopting EOS_GAMMA
//                   --> By linking to the function pointers "EoS_InitCPU/GPUFunc_Ptr"
//                2. Must obtain the CPU and GPU function pointers by separate routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_InitFunc_IdealGas( CPU_FuncPtr, GPU_FuncPtr );
//
//                3. Call-by-reference
//
// Parameter   :  EoS_DensEint2Pres_CPU/GPUPtr : CPU/GPU function pointers to be set
//                EoS_DensPres2Eint_CPU/GPUPtr : ...
//                EoS_DensPres2CSqr_CPU/GPUPtr : ...
//
// Return      :  EoS_DensEint2Pres_CPU, EoS_DensPres2Eint_CPU/GPUPtr, EoS_DensPres2CSqr_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void EoS_InitGPUFunc_IdealGas( EoS_DE2P_t &EoS_DensEint2Pres_GPUPtr,
                               EoS_DP2E_t &EoS_DensPres2Eint_GPUPtr,
                               EoS_DP2C_t &EoS_DensPres2CSqr_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Pres_GPUPtr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2Eint_GPUPtr, EoS_DensPres2Eint_Ptr, sizeof(EoS_DP2E_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_InitCPUFunc_IdealGas( EoS_DE2P_t &EoS_DensEint2Pres_CPUPtr,
                               EoS_DP2E_t &EoS_DensPres2Eint_CPUPtr,
                               EoS_DP2C_t &EoS_DensPres2CSqr_CPUPtr )
{
   EoS_DensEint2Pres_CPUPtr = EoS_DensEint2Pres_Ptr;
   EoS_DensPres2Eint_CPUPtr = EoS_DensPres2Eint_Ptr;
   EoS_DensPres2CSqr_CPUPtr = EoS_DensPres2CSqr_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



// =====================================
// III. Set the auxiliary array
// =====================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_InitAuxArray_IdealGas
// Description :  Set the auxiliary array EoS_AuxArray[] for an ideal-gas EoS
//
// Note        :  1. Invoked by EoS_Init() when adopting EOS_GAMMA
//                   --> By linking to the function pointers "EoS_InitAuxArray_Ptr"
//                2. AuxArray[] has the size of EOS_NAUX_MAX defined in Macro.h (default = 10)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Do not change the order of AuxArray[]
//                   --> For example, the dual-energy routines assume AuxArray[0]=GAMMA
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_InitAuxArray_IdealGas( double AuxArray[] )
{

   AuxArray[0] = GAMMA;
   AuxArray[1] = GAMMA - 1.0;
   AuxArray[2] = 1.0 / ( GAMMA - 1.0 );

} // FUNCTION : EoS_InitAuxArray_IdealGas
#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO )
