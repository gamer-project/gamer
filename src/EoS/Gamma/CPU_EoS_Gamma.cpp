#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#include "CUFLU_Shared_FluUtility.cu"
#endif

#if ( MODEL == HYDRO )



/********************************************************
1. Ideal gas EoS with a constant adiabatic index (gamma)

2. This file is shared by both CPU and GPU

   GPU_EoS_Gamma.cu -> CPU_EoS_Gamma.cpp

3. Three steps are required to implement an EoS

   I.   Set an EoS auxiliary array
   II.  Implement EoS conversion functions
   III. Set EoS initialization functions
********************************************************/



// =============================================
// I. Set an EoS auxiliary array
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_SetAuxArray_Gamma
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
//                   AuxArray_Flt[0] = gamma
//                   AuxArray_Flt[1] = gamma-1
//                   AuxArray_Flt[2] = 1/(gamma-1)
//                   AuxArray_Flt[3] = 1/gamma
//
// Note        :  1. Invoked by EoS_Init_Gamma()
//                2. AuxArray_Flt/Int[] have the size of EOS_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Do not change the order of AuxArray_Flt[]
//                   --> For example, the dual-energy routines assume AuxArray_Flt[0]=GAMMA
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_Gamma( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[0] = GAMMA;
   AuxArray_Flt[1] = GAMMA - 1.0;
   AuxArray_Flt[2] = 1.0 / ( GAMMA - 1.0 );
   AuxArray_Flt[3] = 1.0 / GAMMA;

} // FUNCTION : EoS_SetAuxArray_Gamma
#endif // #ifndef __CUDACC__



// =============================================
// II. Implement EoS conversion functions
//     (1) EoS_DensEint2Pres_*
//     (2) EoS_DensPres2Eint_*
//     (3) EoS_DensPres2CSqr_*
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_Gamma
// Description :  Convert gas mass density and internal energy density to gas pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_Gamma() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Pres_Gamma( const real Dens, const real Eint, const real Passive[], const double AuxArray_Flt[],
                                     const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   if ( Hydro_CheckNegative(Dens) )
      printf( "ERROR : invalid input density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Dens, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(Eint) )
      printf( "ERROR : invalid input internal energy (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Eint, __FILE__, __LINE__, __FUNCTION__ );
#  endif // GAMER_DEBUG


   const real Gamma_m1 = (real)AuxArray_Flt[1];
   real Pres;

   Pres = Eint * Gamma_m1;

   return Pres;

} // FUNCTION : EoS_DensEint2Pres_Gamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2Eint_Gamma
// Description :  Convert gas mass density and pressure to gas internal energy density
//
// Note        :  1. See EoS_DensEint2Pres_Gamma()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas internal energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2Eint_Gamma( const real Dens, const real Pres, const real Passive[], const double AuxArray_Flt[],
                                     const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   if ( Hydro_CheckNegative(Dens) )
      printf( "ERROR : invalid input density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Dens, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(Pres) )
      printf( "ERROR : invalid input pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Pres, __FILE__, __LINE__, __FUNCTION__ );
#  endif // GAMER_DEBUG


   const real _Gamma_m1 = (real)AuxArray_Flt[2];
   real Eint;

   Eint = Pres * _Gamma_m1;

   return Eint;

} // FUNCTION : EoS_DensPres2Eint_Gamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_Gamma
// Description :  Convert gas mass density and pressure to sound speed squared
//
// Note        :  1. See EoS_DensEint2Pres_Gamma()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Sound speed square
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_Gamma( const real Dens, const real Pres, const real Passive[], const double AuxArray_Flt[],
                                     const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   if ( Hydro_CheckNegative(Dens) )
      printf( "ERROR : invalid input density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Dens, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(Pres) )
      printf( "ERROR : invalid input pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Pres, __FILE__, __LINE__, __FUNCTION__ );
#  endif // GAMER_DEBUG


   const real Gamma = (real)AuxArray_Flt[0];
   real Cs2;

   Cs2 = Gamma * Pres / Dens;

   return Cs2;

} // FUNCTION : EoS_DensPres2CSqr_Gamma



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_Gamma;
FUNC_SPACE EoS_DP2E_t EoS_DensPres2Eint_Ptr = EoS_DensPres2Eint_Gamma;
FUNC_SPACE EoS_DP2C_t EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_Gamma;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_SetCPU/GPUFunc_Gamma
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_Gamma()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_SetFunc_Gamma( CPU_FuncPtr, GPU_FuncPtr );
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
void EoS_SetGPUFunc_Gamma( EoS_DE2P_t &EoS_DensEint2Pres_GPUPtr,
                           EoS_DP2E_t &EoS_DensPres2Eint_GPUPtr,
                           EoS_DP2C_t &EoS_DensPres2CSqr_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Pres_GPUPtr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2Eint_GPUPtr, EoS_DensPres2Eint_Ptr, sizeof(EoS_DP2E_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_Gamma( EoS_DE2P_t &EoS_DensEint2Pres_CPUPtr,
                           EoS_DP2E_t &EoS_DensPres2Eint_CPUPtr,
                           EoS_DP2C_t &EoS_DensPres2CSqr_CPUPtr )
{
   EoS_DensEint2Pres_CPUPtr = EoS_DensEint2Pres_Ptr;
   EoS_DensPres2Eint_CPUPtr = EoS_DensPres2Eint_Ptr;
   EoS_DensPres2CSqr_CPUPtr = EoS_DensPres2CSqr_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void EoS_SetAuxArray_Gamma( double [], int [] );
void EoS_SetCPUFunc_Gamma( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t & );
#ifdef GPU
void EoS_SetGPUFunc_Gamma( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_Gamma
// Description :  Initialize EoS
//
// Note        :  1. Set an auxiliary array by invoking EoS_SetAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU EoS routines by invoking EoS_SetCPU/GPUFunc_*()
//                3. Invoked by EoS_Init()
//                   --> Enable it by linking to the function pointer "EoS_Init_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void EoS_Init_Gamma()
{

   EoS_SetAuxArray_Gamma( EoS_AuxArray_Flt, EoS_AuxArray_Int );
   EoS_SetCPUFunc_Gamma( EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr, EoS_DensPres2CSqr_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_Gamma( EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr, EoS_DensPres2CSqr_GPUPtr );
#  endif

} // FUNCTION : EoS_Init_Gamma

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO )
