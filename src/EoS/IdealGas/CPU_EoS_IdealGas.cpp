#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#if ( MODEL == HYDRO )



/*****************************************
This file is shared by both CPU and GPU

    CPU_EoS_IdealGas.cpp
    GPU_EoS_IdealGas.cu

*****************************************/



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_IdealGas
// Description :  Convert gas mass density and internal energy density to gas pressure using an ideal-gas EoS
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. Auxiliary array UserArray[] is set by EoS_IdealGas_InitAuxArray_DensEint2Pres(), where
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

   Pres = Eint*Gamma_m1;

   return Pres;

} // FUNCTION : EoS_DensEint2Pres_IdealGas



// =================================
// get the CPU/GPU function pointers
// =================================

#ifdef __CUDACC__
__device__
#endif
static EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_IdealGas;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUEoS_DensEint2Pres_IdealGas
// Description :  Return the function pointers to the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init() when adopting EOS_IDEALGAS
//                   --> By linking to the function pointers "SetCPU/GPUEoS_DensEint2Pres_Ptr"
//
// Parameter   :  CPU/GPUEoS_DensEint2Pres_Ptr (call-by-reference)
//
// Return      :  CPU/GPUEoS_DensEint2Pres_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUEoS_DensEint2Pres_IdealGas( EoS_DE2P_t &GPUEoS_DensEint2Pres_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUEoS_DensEint2Pres_Ptr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUEoS_DensEint2Pres_IdealGas( EoS_DE2P_t &CPUEoS_DensEint2Pres_Ptr )
{
   CPUEoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_InitAuxArray_IdealGas
// Description :  Set the auxiliary array EoS_AuxArray[] for an ideal-gas EoS
//
// Note        :  1. Invoked by EoS_Init() when adopting EOS_IDEALGAS
//                   --> By linking to the function pointers "EoS_InitAuxArray_Ptr"
//                2. AuxArray[] has the size of EOS_NAUX_MAX defined in CUFLU.h (default = 10)
//                3. Add "#ifndef __CUDACC__" since it is only useful in the CPU code
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
