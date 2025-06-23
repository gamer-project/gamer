#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#include "CUFLU_Shared_FluUtility.cu"
#endif

#if ( MODEL == HYDRO )

extern double rho_AD_SinkParTest;


/********************************************************
1. Template of a user-defined EoS (EOS_USER)

2. This file is shared by both CPU and GPU

   GPU_EoS_Barotropic_SinkParTest.cu -> CPU_EoS_Barotropic_SinkParTest.cpp

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
// Function    :  EoS_SetAuxArray_Barotropic_SinkParTest
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by EoS_Init_Barotropic_SinkParTest()
//                2. AuxArray_Flt/Int[] have the size of EOS_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Physical constants such as Const_amu/Const_kB should be set to unity when disabling OPT__UNIT
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_Barotropic_SinkParTest( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[0] = GAMMA;
   AuxArray_Flt[1] = GAMMA - 1.0;
   AuxArray_Flt[2] = 1.0 / ( GAMMA - 1.0 );
   AuxArray_Flt[3] = 1.0 / GAMMA;
   AuxArray_Flt[4] = ( OPT__UNIT ) ? MOLECULAR_WEIGHT * MU_NORM / Const_kB * (UNIT_E/UNIT_M)
                                   : MOLECULAR_WEIGHT;
   AuxArray_Flt[5] = 1.0 / AuxArray_Flt[4];
   AuxArray_Flt[6] = ISO_TEMP;
   AuxArray_Flt[7] = rho_AD_SinkParTest;

   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "   Background temperature       = %13.7e K\n"    ,    ISO_TEMP      );
      Aux_Message( stdout, "   Adiabatic density            = %13.7e g/cm3\n",    rho_AD_SinkParTest*UNIT_D );
   }

} // FUNCTION : EoS_SetAuxArray_Barotropic_SinkParTest
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
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_Barotropic_SinkParTest
// Description :  Convert gas mass density and internal energy density to gas pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_Barotropic_SinkParTest() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Pres_Barotropic_SinkParTest( const real Dens, const real Eint, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  if ( NCOMP_PASSIVE > 0 )
   if ( Passive == NULL )  printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
#  endif
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical( UNPHY_MODE_SING, &Dens, "input density",
                       TINY_NUMBER, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                       ERROR_INFO, UNPHY_VERBOSE );
// note that some EoS may support Eint<0
   // Hydro_IsUnphysical( UNPHY_MODE_SING, &Eint, "input internal energy",
   //                     (real)0.0,   HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
   //                     ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real Gamma_m1  = (real)AuxArray_Flt[1];
   const real _m_kB     = (real)AuxArray_Flt[5];
   const real T0        = (real)AuxArray_Flt[6];
   const real rho_AD_SinkParTest    = (real)AuxArray_Flt[7];
   real Temp, Cs2, Pres;

   Temp = T0*( 1+ POW( Dens / rho_AD_SinkParTest, Gamma_m1 ) );
   Cs2  = Temp*_m_kB;
   Pres = Cs2*Dens;


// check
#  ifdef GAMER_DEBUG
   if (  Hydro_IsUnphysical( UNPHY_MODE_SING, &Pres, "output pressure",
                             (real)0.0, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                             ERROR_INFO, UNPHY_VERBOSE )  )
   {
      printf( "Dens=%13.7e, Eint=%13.7e\n", Dens, Eint );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Pres;

} // FUNCTION : EoS_DensEint2Pres_Barotropic_SinkParTest



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2Eint_Barotropic_SinkParTest
// Description :  Convert gas mass density and pressure to gas internal energy density
//
// Note        :  1. See EoS_DensEint2Pres_Barotropic_SinkParTest()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas internal energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2Eint_Barotropic_SinkParTest( const real Dens, const real Pres, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  if ( NCOMP_PASSIVE > 0 )
   if ( Passive == NULL )  printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
#  endif
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical( UNPHY_MODE_SING, &Dens, "input density",
                       TINY_NUMBER, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                       ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical( UNPHY_MODE_SING, &Pres, "input pressure",
                       (real)0.0,   HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                       ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real _Gamma_m1 = (real)AuxArray_Flt[2];
   real Eint;

   Eint = Pres * _Gamma_m1;


// check
#  ifdef GAMER_DEBUG
// note that some EoS may support Eint<0
   if (  Hydro_IsUnphysical( UNPHY_MODE_SING, &Eint, "output internal energy",
                             (real)0.0, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                             ERROR_INFO, UNPHY_VERBOSE )  )
   {
      printf( "Dens=%13.7e, Pres=%13.7e\n", Dens, Pres );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Eint;

} // FUNCTION : EoS_DensPres2Eint_Barotropic_SinkParTest



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_Barotropic_SinkParTest
// Description :  Convert gas mass density and pressure to sound speed squared
//
// Note        :  1. See EoS_DensEint2Pres_Barotropic_SinkParTest()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Sound speed squared
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_Barotropic_SinkParTest( const real Dens, const real Pres, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  if ( NCOMP_PASSIVE > 0 )
   if ( Passive == NULL )  printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
#  endif
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical( UNPHY_MODE_SING, &Dens, "input density",
                       TINY_NUMBER, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                       ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical( UNPHY_MODE_SING, &Pres, "input pressure",
                       (real)0.0,   HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                       ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG

   const real Gamma_m1  = (real)AuxArray_Flt[1];
   const real _m_kB     = (real)AuxArray_Flt[5];
   const real T0        = (real)AuxArray_Flt[6];
   const real rho_AD_SinkParTest    = (real)AuxArray_Flt[7];
   real Temp, Cs2;

   Temp = T0*( 1+ POW( Dens / rho_AD_SinkParTest, Gamma_m1 ) );
   Cs2  = Temp*_m_kB;

// check
#  ifdef GAMER_DEBUG
   if (  Hydro_IsUnphysical( UNPHY_MODE_SING, &Cs2, "output sound speed squared",
                             (real)0.0, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                             ERROR_INFO, UNPHY_VERBOSE )  )
   {
      printf( "Dens=%13.7e, Pres=%13.7e\n", Dens, Pres );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Cs2;

} // FUNCTION : EoS_DensPres2CSqr_Barotropic_SinkParTest



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Temp_Barotropic_SinkParTest
// Description :  Convert gas mass density and internal energy density to gas temperature
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_Barotropic_SinkParTest() for the values stored in AuxArray_Flt/Int[]
//                3. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas temperature in kelvin
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Temp_Barotropic_SinkParTest( const real Dens, const real Eint, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical( UNPHY_MODE_SING, &Dens, "input density",
                       TINY_NUMBER, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                       ERROR_INFO, UNPHY_VERBOSE );
// note that some EoS may support Eint<0
   // Hydro_IsUnphysical( UNPHY_MODE_SING, &Eint, "input internal energy",
   //                     (real)0.0,   HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
   //                     ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG

   const real Gamma_m1  = (real)AuxArray_Flt[1];
   const real T0        = (real)AuxArray_Flt[6];
   const real rho_AD_SinkParTest    = (real)AuxArray_Flt[7];
   real Temp;

   Temp = T0*( 1+ POW( Dens / rho_AD_SinkParTest, Gamma_m1 ) );


// check
#  ifdef GAMER_DEBUG
   if (  Hydro_IsUnphysical( UNPHY_MODE_SING, &Temp, "output temperature",
                             (real)0.0, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                             ERROR_INFO, UNPHY_VERBOSE )  )
   {
      printf( "Dens=%13.7e, Eint=%13.7e\n", Dens, Eint );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Temp;

} // FUNCTION : EoS_DensEint2Temp_Barotropic_SinkParTest



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensTemp2Pres_Barotropic_SinkParTest
// Description :  Convert gas mass density and temperature to gas pressure
//
// Note        :  1. See EoS_SetAuxArray_Barotropic_SinkParTest() for the values stored in AuxArray_Flt/Int[]
//                2. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Temp       : Gas temperature in kelvin
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensTemp2Pres_Barotropic_SinkParTest( const real Dens, const real Temp, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical( UNPHY_MODE_SING, &Dens, "input density",
                       TINY_NUMBER, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                       ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical( UNPHY_MODE_SING, &Temp, "input temperature",
                       (real)0.0,   HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                       ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG

   const real _m_kB     = (real)AuxArray_Flt[5];
   real Cs2, Pres;

   Cs2  = Temp*_m_kB;
   Pres = Cs2*Dens;


// check
#  ifdef GAMER_DEBUG
   if (  Hydro_IsUnphysical( UNPHY_MODE_SING, &Pres, "output pressure",
                             (real)0.0, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                             ERROR_INFO, UNPHY_VERBOSE )  )
   {
      printf( "Dens=%13.7e, Temp=%13.7e\n", Dens, Temp );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Pres;

} // FUNCTION : EoS_DensTemp2Pres_Barotropic_SinkParTest



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Entr_Barotropic_SinkParTest
// Description :  Convert gas mass density and internal energy density to gas entropy
//
// Note        :  1. See EoS_SetAuxArray_Barotropic_SinkParTest() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas entropy
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Entr_Barotropic_SinkParTest( const real Dens, const real Eint, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical( UNPHY_MODE_SING, &Dens, "input density",
                       TINY_NUMBER, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
                       ERROR_INFO, UNPHY_VERBOSE );
// note that some EoS may support Eint<0
   // Hydro_IsUnphysical( UNPHY_MODE_SING, &Eint, "input internal energy",
   //                     (real)0.0,   HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
   //                     ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   // real Entr = -1.0;

   /*
   Entr = ...;
   */


// check
// #  ifdef GAMER_DEBUG
//    if (  Hydro_IsUnphysical( UNPHY_MODE_SING, &Entr, "output entropy",
//                              (real)0.0, HUGE_NUMBER, NULL_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
//                              ERROR_INFO, UNPHY_VERBOSE )  )
//    {
//       printf( "Dens=%13.7e, Eint=%13.7e\n", Dens, Eint );
// #     if ( NCOMP_PASSIVE > 0 )
//       printf( "Passive scalars:" );
//       for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
//       printf( "\n" );
// #     endif
//    }
// #  endif // GAMER_DEBUG


    return NULL_REAL;

} // FUNCTION : EoS_DensEint2Entr_Barotropic_SinkParTest



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_General_Barotropic_SinkParTest
// Description :  General EoS converter: In_*[] -> Out[]
//
// Note        :  1. See EoS_DensEint2Pres_Barotropic_SinkParTest()
//                2. In_*[] and Out[] must NOT overlap
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
static void EoS_General_Barotropic_SinkParTest( const int Mode, real Out[], const real In_Flt[], const int In_Int[],
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( Out          == NULL )   printf( "ERROR : Out == NULL in %s !!\n", __FUNCTION__ );
   if ( In_Flt       == NULL )   printf( "ERROR : In_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );
#  endif // GAMER_DEBUG


   /*
   if ( Mode == ... )   Out[...] = ...;
   */

} // FUNCTION : EoS_General_Barotropic_SinkParTest



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_Barotropic_SinkParTest;
FUNC_SPACE EoS_DP2E_t EoS_DensPres2Eint_Ptr = EoS_DensPres2Eint_Barotropic_SinkParTest;
FUNC_SPACE EoS_DP2C_t EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_Barotropic_SinkParTest;
FUNC_SPACE EoS_DE2T_t EoS_DensEint2Temp_Ptr = EoS_DensEint2Temp_Barotropic_SinkParTest;
FUNC_SPACE EoS_DT2P_t EoS_DensTemp2Pres_Ptr = EoS_DensTemp2Pres_Barotropic_SinkParTest;
FUNC_SPACE EoS_DE2S_t EoS_DensEint2Entr_Ptr = EoS_DensEint2Entr_Barotropic_SinkParTest;
FUNC_SPACE EoS_GENE_t EoS_General_Ptr       = EoS_General_Barotropic_SinkParTest;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_SetCPU/GPUFunc_Barotropic_SinkParTest
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_Barotropic_SinkParTest()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_SetFunc_Barotropic_SinkParTest( CPU_FuncPtr, GPU_FuncPtr );
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
//
// Return      :  EoS_DensEint2Pres_CPU/GPUPtr, EoS_DensPres2Eint_CPU/GPUPtr,
//                EoS_DensPres2CSqr_CPU/GPUPtr, EoS_DensEint2Temp_CPU/GPUPtr,
//                EoS_DensTemp2Pres_CPU/GPUPtr, EoS_DensEint2Entr_CPU/GPUPtr,
//                EoS_General_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void EoS_SetGPUFunc_Barotropic_SinkParTest( EoS_DE2P_t &EoS_DensEint2Pres_GPUPtr,
                                   EoS_DP2E_t &EoS_DensPres2Eint_GPUPtr,
                                   EoS_DP2C_t &EoS_DensPres2CSqr_GPUPtr,
                                   EoS_DE2T_t &EoS_DensEint2Temp_GPUPtr,
                                   EoS_DT2P_t &EoS_DensTemp2Pres_GPUPtr,
                                   EoS_DE2S_t &EoS_DensEint2Entr_GPUPtr,
                                   EoS_GENE_t &EoS_General_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Pres_GPUPtr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2Eint_GPUPtr, EoS_DensPres2Eint_Ptr, sizeof(EoS_DP2E_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Temp_GPUPtr, EoS_DensEint2Temp_Ptr, sizeof(EoS_DE2T_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensTemp2Pres_GPUPtr, EoS_DensTemp2Pres_Ptr, sizeof(EoS_DT2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Entr_GPUPtr, EoS_DensEint2Entr_Ptr, sizeof(EoS_DE2S_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_General_GPUPtr,       EoS_General_Ptr,       sizeof(EoS_GENE_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_Barotropic_SinkParTest( EoS_DE2P_t &EoS_DensEint2Pres_CPUPtr,
                                   EoS_DP2E_t &EoS_DensPres2Eint_CPUPtr,
                                   EoS_DP2C_t &EoS_DensPres2CSqr_CPUPtr,
                                   EoS_DE2T_t &EoS_DensEint2Temp_CPUPtr,
                                   EoS_DT2P_t &EoS_DensTemp2Pres_CPUPtr,
                                   EoS_DE2S_t &EoS_DensEint2Entr_CPUPtr,
                                   EoS_GENE_t &EoS_General_CPUPtr )
{
   EoS_DensEint2Pres_CPUPtr = EoS_DensEint2Pres_Ptr;
   EoS_DensPres2Eint_CPUPtr = EoS_DensPres2Eint_Ptr;
   EoS_DensPres2CSqr_CPUPtr = EoS_DensPres2CSqr_Ptr;
   EoS_DensEint2Temp_CPUPtr = EoS_DensEint2Temp_Ptr;
   EoS_DensTemp2Pres_CPUPtr = EoS_DensTemp2Pres_Ptr;
   EoS_DensEint2Entr_CPUPtr = EoS_DensEint2Entr_Ptr;
   EoS_General_CPUPtr       = EoS_General_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void EoS_SetAuxArray_Barotropic_SinkParTest( double [], int [] );
void EoS_SetCPUFunc_Barotropic_SinkParTest( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_DE2S_t &, EoS_GENE_t & );
#ifdef GPU
void EoS_SetGPUFunc_Barotropic_SinkParTest( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_DE2S_t &, EoS_GENE_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_Barotropic_SinkParTest
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
void EoS_Init_Barotropic_SinkParTest()
{

   EoS_SetAuxArray_Barotropic_SinkParTest( EoS_AuxArray_Flt, EoS_AuxArray_Int );
   EoS_SetCPUFunc_Barotropic_SinkParTest( EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                                 EoS_DensPres2CSqr_CPUPtr, EoS_DensEint2Temp_CPUPtr,
                                 EoS_DensTemp2Pres_CPUPtr, EoS_DensEint2Entr_CPUPtr,
                                 EoS_General_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_Barotropic_SinkParTest( EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr,
                                 EoS_DensPres2CSqr_GPUPtr, EoS_DensEint2Temp_GPUPtr,
                                 EoS_DensTemp2Pres_GPUPtr, EoS_DensEint2Entr_GPUPtr,
                                 EoS_General_GPUPtr );
#  endif

} // FUNCTION : EoS_Init_Barotropic_SinkParTest

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO )
