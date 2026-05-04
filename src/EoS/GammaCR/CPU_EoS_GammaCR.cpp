#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#include "CUFLU_Shared_FluUtility.cu"
#endif

#if ( MODEL == HYDRO )

#ifdef COSMIC_RAY

#ifdef __CUDACC__
__device__ static real EoS_CREint2CRPres_GammaCR( const real E_CR,
                                                  const double AuxArray_Flt[], const int AuxArray_Int[],
                                                  const real *const Table[EOS_NTABLE_MAX] );
#else // #ifdef __CUDACC__
static real EoS_CREint2CRPres_GammaCR( const real E_CR,
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX] );
#endif // #ifdef __CUDACC__ ... else ...


/********************************************************
1. Ideal gas EoS with cosmic-ray component (GAMMA_CR)

2. This file is shared by both CPU and GPU

   GPU_EoS_GammaCR.cu -> CPU_EoS_GammaCR.cpp

3. Three steps are required to implement an EoS

   I.   Set EoS auxiliary arrays
   II.  Implement EoS conversion functions
   III. Set EoS initialization functions

4. All EoS conversion functions must be thread-safe and
   not use any global variable

5. When an EoS conversion function fails, it is recommended
   to return NAN in order to trigger auto-correction such as
   "OPT__1ST_FLUX_CORR" and "AUTO_REDUCE_DT"
********************************************************/



// =============================================
// I. Set EoS auxiliary arrays
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_SetAuxArray_GammaCR
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
//                   AuxArray_Flt[0] = gamma_gas
//                   AuxArray_Flt[1] = gamma_gas - 1
//                   AuxArray_Flt[2] = 1/(gamma_gas - 1)
//                   AuxArray_Flt[3] = 1/gamma_gas
//                   AuxArray_Flt[4] = gamma_cr
//                   AuxArray_Flt[5] = gamma_cr - 1
//                   AuxArray_Flt[6] = minimum pressure
//                   AuxArray_Flt[7] = (mean molecular weight)*(atomic mass unit)/(Boltzmann constant)*(UNIT_E/UNIT_M)
//                   AuxArray_Flt[8] = 1/AuxArray_Flt[7]
//
// Note        :  1. Invoked by EoS_Init_GammaCR()
//                2. AuxArray_Flt/Int[] have the size of EOS_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Physical constants such as MU_NORM/Const_kB should be set to unity when disabling OPT__UNIT
//                5. Do not change the order of AuxArray_Flt[]
//                   --> For example, the dual-energy routines assume AuxArray_Flt[0]=GAMMA
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_GammaCR( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[0] = GAMMA;
   AuxArray_Flt[1] = GAMMA - 1.0;
   AuxArray_Flt[2] = 1.0 / ( GAMMA - 1.0);
   AuxArray_Flt[3] = 1.0 / GAMMA;
   AuxArray_Flt[4] = GAMMA_CR;
   AuxArray_Flt[5] = GAMMA_CR - 1.0;
   AuxArray_Flt[6] = MIN_PRES;
   AuxArray_Flt[7] = ( OPT__UNIT ) ? MOLECULAR_WEIGHT * MU_NORM / Const_kB * (UNIT_E/UNIT_M)
                                   : MOLECULAR_WEIGHT;
   AuxArray_Flt[8] = 1.0 / AuxArray_Flt[7];

} // FUNCTION : EoS_SetAuxArray_GammaCR
#endif // #ifndef __CUDACC__



// =============================================
// II. Implement EoS conversion functions
//     (1) EoS_DensEint2Pres_*
//     (2) EoS_DensPres2Eint_*
//     (3) EoS_DensPres2CSqr_*
//     (4) EoS_DensEint2Temp_* [OPTIONAL]
//     (5) EoS_DensTemp2Pres_* [OPTIONAL]
//     (6) EoS_DensEint2Entr_* [OPTIONAL]
//     (7) EoS_General_*       [OPTIONAL]
//     (8) EoS_CREint2CRPres_*
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_GammaCR
// Description :  Convert gas mass density and total internal energy density (gas + cosmic ray) to total pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_GammaCR() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Total internal energy density (gas + cosmic ray)
//                Passive    : Passive scalars (Passive[CRAY-NCOMP_FLUID] gives the cosmic-ray energy density)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Total pressure (gas + cosmic ray)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Pres_GammaCR( const real Dens, const real Eint, const real Passive[],
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX])
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density",         TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( Eint, "input internal energy", (real)0.0,   HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real Gamma_m1  = (real)AuxArray_Flt[1];
   const real small_val = (real)AuxArray_Flt[6];
   const real E_CR      = Passive[ CRAY-NCOMP_FLUID ];
   const real Pres_CR   = EoS_CREint2CRPres_GammaCR( E_CR, AuxArray_Flt, AuxArray_Int, Table );
   real Pres;

   Pres = Gamma_m1*FMAX( Eint-E_CR, small_val ) + Pres_CR;

   return Pres;

} // FUNCTION : EoS_DensEint2Pres_GammaCR



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2Eint_GammaCR
// Description :  Convert gas mass density and total pressure to total internal energy density (gas + cosmic ray)
//
// Note        :  1. See EoS_DensEint2Pres_GammaCR()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Total pressure (gas + cosmic ray)
//                Passive    : Passive scalars (Passive[CRAY-NCOMP_FLUID] gives the cosmic-ray energy density)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Total internal energy density (gas + cosmic ray)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2Eint_GammaCR( const real Dens, const real Pres, const real Passive[],
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density",  TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( Pres, "input pressure", (real)0.0,   HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real Gamma_m1_inv = (real)AuxArray_Flt[2];
   const real small_val    = (real)AuxArray_Flt[6];
   const real E_CR         = Passive[ CRAY-NCOMP_FLUID ];
   const real Pres_CR      = EoS_CREint2CRPres_GammaCR( E_CR, AuxArray_Flt, AuxArray_Int, Table );
   real Eint;

   Eint = FMAX( Pres-Pres_CR, small_val )*Gamma_m1_inv + E_CR;

   return Eint;

} // FUNCTION : EoS_DensPres2Eint_GammaCR



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_GammaCR
// Description :  Convert gas mass density and total pressure to effective sound speed squared
//
// Note        :  1. See EoS_DensEint2Pres_GammaCR()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Total pressure (gas + cosmic ray)
//                Passive    : Passive scalars (Passive[CRAY-NCOMP_FLUID] gives the cosmic-ray energy density)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Effective sound speed squared
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_GammaCR( const real Dens, const real Pres, const real Passive[],
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density",  TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( Pres, "input pressure", (real)0.0,   HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real Gamma     = (real)AuxArray_Flt[0];
   const real GammaCR   = (real)AuxArray_Flt[4];
   const real small_val = (real)AuxArray_Flt[6];
   const real E_CR      = Passive[ CRAY-NCOMP_FLUID ];
   const real Pres_CR   = EoS_CREint2CRPres_GammaCR( E_CR, AuxArray_Flt, AuxArray_Int, Table );
   real Cs2;

   Cs2 = (  GammaCR*Pres_CR + Gamma*FMAX( Pres-Pres_CR, small_val )  ) / Dens;

   return Cs2;

} // FUNCTION : EoS_DensPres2CSqr_GammaCR



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Temp_GammaCR
// Description :  Convert gas mass density and total internal energy density (gas + cosmic ray) to gas temperature
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_GammaCR() for the values stored in AuxArray_Flt/Int[]
//                3. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Total internal energy density (gas + cosmic ray)
//                Passive    : Passive scalars (Passive[CRAY-NCOMP_FLUID] gives the cosmic-ray energy density)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas temperature in kelvin
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Temp_GammaCR( const real Dens, const real Eint, const real Passive[],
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density",         TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( Eint, "input internal energy", (real)0.0,   HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG

   const real Gamma_m1  = (real)AuxArray_Flt[1];
   const real small_val = (real)AuxArray_Flt[6];
   const real m_kB      = (real)AuxArray_Flt[7];
   const real E_CR      = Passive[ CRAY-NCOMP_FLUID ];
   real Pres_Gas, Temp;

   Pres_Gas = FMAX( Eint-E_CR, small_val )*Gamma_m1;
   Temp     = m_kB*Pres_Gas/Dens;

   return Temp;

} // FUNCTION : EoS_DensEint2Temp_GammaCR


//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensTemp2Pres_GammaCR
// Description :  Convert gas mass density and temperature to total pressure
//
// Note        :  1. See EoS_SetAuxArray_GammaCR() for the values stored in AuxArray_Flt/Int[]
//                2. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Temp       : Gas temperature in kelvin
//                Passive    : Passive scalars (Passive[CRAY-NCOMP_FLUID] gives the cosmic-ray energy density)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Total pressure (gas + cosmic ray)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensTemp2Pres_GammaCR( const real Dens, const real Temp, const real Passive[],
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density",     TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( Temp, "input temperature", (real)0.0,   HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real _m_kB   = (real)AuxArray_Flt[8];
   const real E_CR    = Passive[ CRAY-NCOMP_FLUID ];
   const real Pres_CR = EoS_CREint2CRPres_GammaCR( E_CR, AuxArray_Flt, AuxArray_Int, Table );
   real Pres;

   Pres = Temp*Dens*_m_kB + Pres_CR;

   return Pres;

} // FUNCTION : EoS_DensTemp2Pres_GammaCR



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Entr_GammaCR
// Description :  Convert gas mass density and total internal energy density (gas + cosmic ray) to gas entropy
//                --> Here entropy is defined as "pressure_gas / density^(Gamma_gas-1)" (i.e., entropy per volume)
//
// Note        :  1. See EoS_SetAuxArray_GammaCR() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Total internal energy density (gas + cosmic ray)
//                Passive    : Passive scalars (Passive[CRAY-NCOMP_FLUID] gives the cosmic-ray energy density)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas entropy
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Entr_GammaCR( const real Dens, const real Eint, const real Passive[],
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density",         TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( Eint, "input internal energy", (real)0.0,   HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real Gamma_m1  = (real)AuxArray_Flt[1];
   const real small_val = (real)AuxArray_Flt[6];
   const real E_CR      = Passive[ CRAY-NCOMP_FLUID ];
   real Pres_Gas, Entr_Gas;

   Pres_Gas = Gamma_m1*FMAX( Eint-E_CR, small_val );
   Entr_Gas = Pres_Gas*POW( Dens, -Gamma_m1 );

   return Entr_Gas;

} // FUNCTION : EoS_DensEint2Entr_GammaCR



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_General_GammaCR
// Description :  General EoS converter: In[] -> Out[]
//
// Note        :  1. See EoS_DensEint2Pres_GammaCR()
//                2. In_*[] and Out[] must NOT overlap
//                3. Useless for this EoS
//
// Parameter   :  Mode       : To support multiple modes in this general converter
//                Out        : Output array
//                In_*       : Input array
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Out[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void EoS_General_GammaCR( const int Mode, real Out[], const real In_Flt[], const int In_Int[],
                                 const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real *const Table[EOS_NTABLE_MAX] )
{

// not used by this EoS

} // FUNCTION : EoS_General_GammaCR



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_CREint2CRPres_GammaCR
// Description :  Convert cosmic-ray energy density to cosmic-ray pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_GammaCR() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  E_CR       : Cosmic-ray energy density
//                AuxArray_* : Auxiliary arrays (see the Note above)
//
// Return      :  Cosmic ray pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_CREint2CRPres_GammaCR( const real E_CR,
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( E_CR < (real)0.0 )
      printf( "ERROR : invalid input cosmic-ray energy density (%13.7e) in %s() !!\n", E_CR, __FUNCTION__ );
#  endif // GAMER_DEBUG


   const real GammaCR_m1 = (real)AuxArray_Flt[5];
   real Pres_CR;

   Pres_CR = GammaCR_m1*E_CR;


// check
#  ifdef GAMER_DEBUG
   if ( Pres_CR < (real)0.0 )
   {
      printf( "ERROR : invalid output cosmic-ray pressure (%13.7e) in %s() !!\n", Pres_CR, __FUNCTION__ );
      printf( "        CRay=%13.7e\n", E_CR );
   }
#  endif // GAMER_DEBUG


   return Pres_CR;

} // FUNCTION : EoS_CREint2CRPres_GammaCR



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_DE2P_t    EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_GammaCR;
FUNC_SPACE EoS_DP2E_t    EoS_DensPres2Eint_Ptr = EoS_DensPres2Eint_GammaCR;
FUNC_SPACE EoS_DP2C_t    EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_GammaCR;
FUNC_SPACE EoS_DE2T_t    EoS_DensEint2Temp_Ptr = EoS_DensEint2Temp_GammaCR;
FUNC_SPACE EoS_DT2P_t    EoS_DensTemp2Pres_Ptr = EoS_DensTemp2Pres_GammaCR;
FUNC_SPACE EoS_DE2S_t    EoS_DensEint2Entr_Ptr = EoS_DensEint2Entr_GammaCR;
FUNC_SPACE EoS_GENE_t    EoS_General_Ptr       = EoS_General_GammaCR;
FUNC_SPACE EoS_CRE2CRP_t EoS_CREint2CRPres_Ptr = EoS_CREint2CRPres_GammaCR;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_SetCPU/GPUFunc_GammaCR
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_GammaCR()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_SetFunc_GammaCR( CPU_FuncPtr, GPU_FuncPtr );
//
//                3. Call-by-reference
//
// Parameter   :  EoS_DensEint2Pres_CPU/GPUPtr : CPU/GPU function pointers to be set
//                EoS_DensPres2Eint_CPU/GPUPtr : ...
//                EoS_DensPres2CSqr_CPU/GPUPtr : ...
//                EoS_DensEint2Temp_CPU/GPUPtr : ...
//                EoS_DensTemp2Pres_CPU/GPUPtr : ...
//                EoS_DensEint2Entr_CPU/GPUPtr : ...
//                EoS_General_CPU/GPUPtr       : ...
//                EoS_CREint2CRPres_CPU/GPUPtr : ...
//
// Return      :  EoS_DensEint2Pres_CPU/GPUPtr, EoS_DensPres2Eint_CPU/GPUPtr,
//                EoS_DensPres2CSqr_CPU/GPUPtr, EoS_DensEint2Temp_CPU/GPUPtr,
//                EoS_DensTemp2Pres_CPU/GPUPtr, EoS_DensEint2Entr_CPU/GPUPtr,
//                EoS_General_CPU/GPUPtr, EoS_CREint2CRPres_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void EoS_SetGPUFunc_GammaCR( EoS_DE2P_t    &EoS_DensEint2Pres_GPUPtr,
                             EoS_DP2E_t    &EoS_DensPres2Eint_GPUPtr,
                             EoS_DP2C_t    &EoS_DensPres2CSqr_GPUPtr,
                             EoS_DE2T_t    &EoS_DensEint2Temp_GPUPtr,
                             EoS_DT2P_t    &EoS_DensTemp2Pres_GPUPtr,
                             EoS_DE2S_t    &EoS_DensEint2Entr_GPUPtr,
                             EoS_GENE_t    &EoS_General_GPUPtr,
                             EoS_CRE2CRP_t &EoS_CREint2CRPres_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Pres_GPUPtr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t   ) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2Eint_GPUPtr, EoS_DensPres2Eint_Ptr, sizeof(EoS_DP2E_t   ) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t   ) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Temp_GPUPtr, EoS_DensEint2Temp_Ptr, sizeof(EoS_DE2T_t   ) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensTemp2Pres_GPUPtr, EoS_DensTemp2Pres_Ptr, sizeof(EoS_DT2P_t   ) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Entr_GPUPtr, EoS_DensEint2Entr_Ptr, sizeof(EoS_DE2S_t   ) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_General_GPUPtr,       EoS_General_Ptr,       sizeof(EoS_GENE_t   ) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_CREint2CRPres_GPUPtr, EoS_CREint2CRPres_Ptr, sizeof(EoS_CRE2CRP_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_GammaCR( EoS_DE2P_t    &EoS_DensEint2Pres_CPUPtr,
                             EoS_DP2E_t    &EoS_DensPres2Eint_CPUPtr,
                             EoS_DP2C_t    &EoS_DensPres2CSqr_CPUPtr,
                             EoS_DE2T_t    &EoS_DensEint2Temp_CPUPtr,
                             EoS_DT2P_t    &EoS_DensTemp2Pres_CPUPtr,
                             EoS_DE2S_t    &EoS_DensEint2Entr_CPUPtr,
                             EoS_GENE_t    &EoS_General_CPUPtr,
                             EoS_CRE2CRP_t &EoS_CREint2CRPres_CPUPtr )
{
   EoS_DensEint2Pres_CPUPtr = EoS_DensEint2Pres_Ptr;
   EoS_DensPres2Eint_CPUPtr = EoS_DensPres2Eint_Ptr;
   EoS_DensPres2CSqr_CPUPtr = EoS_DensPres2CSqr_Ptr;
   EoS_DensEint2Temp_CPUPtr = EoS_DensEint2Temp_Ptr;
   EoS_DensTemp2Pres_CPUPtr = EoS_DensTemp2Pres_Ptr;
   EoS_DensEint2Entr_CPUPtr = EoS_DensEint2Entr_Ptr;
   EoS_General_CPUPtr       = EoS_General_Ptr;
   EoS_CREint2CRPres_CPUPtr = EoS_CREint2CRPres_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void EoS_SetAuxArray_GammaCR( double [] , int []);
void EoS_SetCPUFunc_GammaCR( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_DE2S_t &, EoS_GENE_t &, EoS_CRE2CRP_t & );
#ifdef GPU
void EoS_SetGPUFunc_GammaCR( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_DE2S_t &, EoS_GENE_t &, EoS_CRE2CRP_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_GammaCR
// Description :  Initialize EoS
//
// Note        :  1. Set auxiliary arrays by invoking EoS_SetAuxArray_*()
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
void EoS_Init_GammaCR()
{

   EoS_SetAuxArray_GammaCR( EoS_AuxArray_Flt, EoS_AuxArray_Int );
   EoS_SetCPUFunc_GammaCR( EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                           EoS_DensPres2CSqr_CPUPtr, EoS_DensEint2Temp_CPUPtr,
                           EoS_DensTemp2Pres_CPUPtr, EoS_DensEint2Entr_CPUPtr,
                           EoS_General_CPUPtr, EoS_CREint2CRPres_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_GammaCR( EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr,
                           EoS_DensPres2CSqr_GPUPtr, EoS_DensEint2Temp_GPUPtr,
                           EoS_DensTemp2Pres_GPUPtr, EoS_DensEint2Entr_GPUPtr,
                           EoS_General_GPUPtr, EoS_CREint2CRPres_GPUPtr );
#  endif

} // FUNCTION : EoS_Init_GammaCR

#endif // #ifndef __CUDACC__



#endif // #ifdef COSMIC_RAY

#endif // #if ( MODEL == HYDRO )
