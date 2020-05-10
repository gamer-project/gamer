#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#if ( MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_IdealGas_DensEint2Pres
// Description :  Convert gas mass density and internal energy density to gas pressure using an ideal-gas EoS
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Internal energy density here is per unit volume instead of per unit mass
//                3. Auxiliary array UserArray[] is set by EoS_IdealGas_InitAuxArray_DensEint2Pres(), where
//                   UserArray[0] = gamma - 1
//
// Parameter   :  Dens      : Gas mass density
//                Eint      : Gas internal energy density
//                UserArray : User-provided auxiliary array (see the Note above)
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_IdealGas_DensEint2Pres( const real Dens, const real Eint, const real UserArray[] )
{

   const real Gamma_m1 = UserArray[0];
   real Pres;

   Pres = Eint*Gamma_m1;

   return Pres;

} // FUNCTION : EoS_IdealGas_DensEint2Pres



// =================================
// get the CPU/GPU function pointers
// =================================

#ifdef __CUDACC__
__device__
#endif
static EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_IdealGas_DensEint2Pres;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_PointMass
// Description :  Return the function pointers to the CPU/GPU EoS routines
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
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
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



#endif // #if ( MODEL == HYDRO )
