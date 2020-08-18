#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#if ( MODEL == HYDRO )



/********************************************************
1. Template of a user-defined EoS (EOS_USER)

2. This file is shared by both CPU and GPU

   GPU_EoS_User_Template.cu -> CPU_EoS_User_Template.cpp

3. Three steps are required to implement an EoS

   I.   Implement EoS conversion functions
   II.  Set an EoS auxiliary array
   III. Set EoS initialization functions
********************************************************/



// =============================================
// I. Implement EoS conversion functions
//    (1) EoS_DensEint2Pres_*
//    (2) EoS_DensPres2Eint_*
//    (3) EoS_DensPres2CSqr_*
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_User_Template
// Description :  Convert gas mass density and internal energy density to gas pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_User_Template() for the values stored in AuxArray[]
//
// Parameter   :  Dens     : Gas mass density
//                Eint     : Gas internal energy density
//                Passive  : Passive scalars
//                AuxArray : Auxiliary array (see the Note above)
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Pres_User_Template( const real Dens, const real Eint, const real Passive[], const double AuxArray[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( Passive  == NULL )    printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray == NULL )    printf( "ERROR : AuxArray == NULL in %s !!\n", __FUNCTION__ );
#  endif

   real Pres = -1.0;

   /*
   Pres = ...;
   */

// check
#  ifdef GAMER_DEBUG
   if ( Pres <= 0.0 )
   {
      printf( "ERROR : invalid pressure (%14.7e) in %s (Dens %14.7e, Eint %14.7e) !!\n",
              Pres, __FUNCTION__, Dens, Eint );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "        Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%14.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif

   return Pres;

} // FUNCTION : EoS_DensEint2Pres_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2Eint_User_Template
// Description :  Convert gas mass density and pressure to gas internal energy density
//
// Note        :  1. See EoS_DensEint2Pres_User_Template()
//
// Parameter   :  Dens     : Gas mass density
//                Pres     : Gas pressure
//                Passive  : Passive scalars
//                AuxArray : Auxiliary array (see the Note above)
//
// Return      :  Gas internal energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2Eint_User_Template( const real Dens, const real Pres, const real Passive[], const double AuxArray[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( Passive  == NULL )    printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray == NULL )    printf( "ERROR : AuxArray == NULL in %s !!\n", __FUNCTION__ );
#  endif

   real Eint = -1.0;

   /*
   Eint = ...;
   */

// check
#  ifdef GAMER_DEBUG
   if ( Eint <= 0.0 )
   {
      printf( "ERROR : invalid internal energy density (%14.7e) in %s (Dens %14.7e, Pres %14.7e) !!\n",
              Eint, __FUNCTION__, Dens, Pres );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "        Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%14.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif

   return Eint;

} // FUNCTION : EoS_DensPres2Eint_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_User_Template
// Description :  Convert gas mass density and pressure to sound speed squared
//
// Note        :  1. See EoS_DensEint2Pres_User_Template()
//
// Parameter   :  Dens     : Gas mass density
//                Pres     : Gas pressure
//                Passive  : Passive scalars
//                AuxArray : Auxiliary array (see the Note above)
//
// Return      :  Sound speed square
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_User_Template( const real Dens, const real Pres, const real Passive[], const double AuxArray[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( Passive  == NULL )    printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray == NULL )    printf( "ERROR : AuxArray == NULL in %s !!\n", __FUNCTION__ );
#  endif

   real Cs2 = -1.0;

   /*
   Cs2 = ...;
   */

// check
#  ifdef GAMER_DEBUG
   if ( Cs2 <= 0.0 )
   {
      printf( "ERROR : invalid sound speed squared (%14.7e) in %s (Dens %14.7e, Pres %14.7e) !!\n",
              Cs2, __FUNCTION__, Dens, Pres );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "        Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%14.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif

   return Cs2;

} // FUNCTION : EoS_DensPres2CSqr_User_Template



// =============================================
// II. Set an EoS auxiliary array
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_SetAuxArray_User_Template
// Description :  Set the auxiliary array AuxArray[]
//
// Note        :  1. Invoked by EoS_Init_User_Template()
//                2. AuxArray[] has the size of EOS_NAUX_MAX defined in Macro.h (default = 10)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_User_Template( double AuxArray[] )
{

   /*
   AuxArray[0] = ...;
   AuxArray[1] = ...;
   */

} // FUNCTION : EoS_SetAuxArray_User_Template
#endif // #ifndef __CUDACC__



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_User_Template;
FUNC_SPACE EoS_DP2E_t EoS_DensPres2Eint_Ptr = EoS_DensPres2Eint_User_Template;
FUNC_SPACE EoS_DP2C_t EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_User_Template;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_InitCPU/GPUFunc_User_Template
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_User_Template()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_InitFunc_User_Template( CPU_FuncPtr, GPU_FuncPtr );
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
void EoS_SetGPUFunc_User_Template( EoS_DE2P_t &EoS_DensEint2Pres_GPUPtr,
                                   EoS_DP2E_t &EoS_DensPres2Eint_GPUPtr,
                                   EoS_DP2C_t &EoS_DensPres2CSqr_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Pres_GPUPtr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2Eint_GPUPtr, EoS_DensPres2Eint_Ptr, sizeof(EoS_DP2E_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_User_Template( EoS_DE2P_t &EoS_DensEint2Pres_CPUPtr,
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
void EoS_SetAuxArray_User_Template( double [] );
void EoS_SetCPUFunc_User_Template( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t & );
#ifdef GPU
void EoS_SetGPUFunc_User_Template( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_User_Template
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
void EoS_Init_User_Template()
{

   EoS_SetAuxArray_User_Template( EoS_AuxArray );
   EoS_SetCPUFunc_User_Template( EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr, EoS_DensPres2CSqr_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_User_Template( EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr, EoS_DensPres2CSqr_GPUPtr );
#  endif

} // FUNCTION : EoS_Init_User_Template

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO )
