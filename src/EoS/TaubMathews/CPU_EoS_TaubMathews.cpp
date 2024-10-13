#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#include "CUFLU_Shared_FluUtility.cu"
#endif

#if ( MODEL == HYDRO  &&  defined SRHD )



/********************************************************
1. Approximately relativistic ideal gas EoS (Taub-Mathews EoS)

2. This file is shared by both CPU and GPU

   GPU_EoS_TaubMathews.cu -> CPU_EoS_TaubMathews.cpp

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
// Function    :  EoS_SetAuxArray_TaubMathews
// Description :  Set the auxiliary array AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by EoS_Init_TaubMathews()
//                2. AuxArray_Flt/Int[] have the size of EOS_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Do not change the order of AuxArray_Flt/Int[]
//                5. Physical constants such as Const_amu/Const_kB should be set to unity when disabling OPT__UNIT
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_TaubMathews( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[0] = ( OPT__UNIT ) ? MOLECULAR_WEIGHT * MU_NORM / Const_kB * (UNIT_E/UNIT_M)
                                   : MOLECULAR_WEIGHT;
   AuxArray_Flt[1] = 1.0 / AuxArray_Flt[0];

} // FUNCTION : EoS_SetAuxArray_TaubMathews
#endif // #ifndef __CUDACC__



// =============================================
// II. Implement EoS conversion functions
//     (1) EoS_HTilde2Temp_*
//     (2) EoS_Temp2HTilde_*
//     (3) EoS_DensPres2CSqr_*
//     (4) EoS_GuessHTilde_*
// =============================================



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_GuessHTilde_TaubMathews
// Description :  Guess a reduced enthalpy for Newton-Raphson iteration
//
// Note        :
//
// Parameter   :  Con         : Conserved variables
//                Constant    : The constant on the LHS of Eq. A3 in "Tseng et al. 2021, MNRAS, 504, 3298"
//                AuxArray_*  : Auxiliary arrays (see the Note above)
//                Table       : EoS tables
//
// Return      :  GuessHTilde : Guessed reduced enthalpy
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_GuessHTilde_TaubMathews( const real Con[], real* const Constant, const double AuxArray_Flt[],
                                         const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

   real GuessHTilde, Discrimination;

   const real Msqr   = SQR(Con[1]) + SQR(Con[2]) + SQR(Con[3]);
   const real Dsqr   = SQR(Con[0]);
   const real abc    = (real)1.0 / Dsqr;
   const real E_D    = Con[4] / Con[0];
   const real M_Dsqr = abc * Msqr;
   const real M_D    = SQRT( M_Dsqr );

// note that replacing a^2-b^2 with (a+b)*(a-b) helps alleviate a catastrophic cancellation
   const real X = SQRT( E_D*E_D + (real)2.0*E_D );
   const real Y = X + M_D;
   const real Z = X - M_D;
   *Constant = Y * Z;

#  ifdef GAMER_DEBUG
   if ( *Constant <= (real)TINY_NUMBER )
      printf( "ERROR : f(HTilde) = %14.7e <= %13.7e (D %13.7e, E_D %13.7e, M_D %13.7e) in %s !!\n",
              *Constant, TINY_NUMBER, Con[0], E_D, M_D, __FUNCTION__ );
#  endif

   const real A = (real)437.0*M_Dsqr + (real)117.0;
   const real B = (real)1.0 + M_Dsqr;

// Eq. A7 in "Tseng et al. 2021, MNRAS, 504, 3298"
   Discrimination  = (real)3240000.0*SQR( B );
   Discrimination /= SQR( A );

   if ( *Constant >= Discrimination )
   {
//    Eq. A6 in "Tseng et al. 2021, MNRAS, 504, 3298"
      GuessHTilde  = (real)4./3.*SQRT( *Constant );
   }

   else
   {
//    Eq. A4 in "Tseng et al. 2021, MNRAS, 504, 3298"
      const real C = (real)43.0*M_Dsqr + (real)63.0;
      const real F = (real)75.0*B;
      const real G = (real)125.0*B*C*(*Constant);
      GuessHTilde  = G;
      GuessHTilde /= C*( F + SQRT(G+F*F) );
   }

   return GuessHTilde;

}



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_HTilde2Temp_TaubMathews
// Description :  Convert reduced enthalpy to temperature
//
// Note        :  Eqs. 16, A2 in "Tseng et al. 2021, MNRAS, 504, 3298"
//
// Parameter   :  HTilde     : Reduced enthalpy
//                Temp       : Temperature (kT/mc**2)
//                DiffTemp   : The derivative of a reduced enthalpy with respect to temperature
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Temp, DiffTemp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void EoS_HTilde2Temp_TaubMathews( const real HTilde, real* const Temp, real* const DiffTemp,
                                         const real Passive[], const double AuxArray_Flt[],
                                         const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

  const real HTildeSqr = SQR(HTilde);
  const real Factor0   = (real)2.0*HTildeSqr + (real)4.0*HTilde;
  const real Factor1   = SQRT( (real)9.0*HTildeSqr + (real)18.0*HTilde + (real)25.0 );
  const real Factor2   = (real)5.0*HTilde + (real)5.0 + Factor1;

  if ( Temp != NULL )
    *Temp = Factor0 / Factor2;

  if ( DiffTemp != NULL )
    *DiffTemp = ( (real)4.0*HTilde + (real)4.0 ) / Factor2
                   -  Factor0*(  ( (real)9.0*HTilde + (real)9.0 ) / Factor1 + (real)5.0  ) / SQR(Factor2);

} // FUNCTION : EoS_HTilde2Temp_TaubMathews



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_Temp2HTilde_TaubMathews
// Description :  Convert temperature to reduced enthalpy
//
// Note        :  Eq. 10 in "Tseng et al. 2021, MNRAS, 504, 3298"
//
// Parameter   :  Temp       : Dimensionless temperature (kT/mc**2)
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Reduced energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_Temp2HTilde_TaubMathews( const real Temp, const real Passive[], const double AuxArray_Flt[],
                                         const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

   const real TempSqr = Temp*Temp;
   const real Factor = (real)2.25*TempSqr;

   real HTilde;   // reduced enthalpy

   HTilde  = Factor;
   HTilde /= (real)1.0 + SQRT( Factor + (real)1.0 );
   HTilde += (real)2.5*Temp;

   return HTilde;

} // FUNCTION : EoS_Temp2HTilde_TaubMathews



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_TaubMathews
// Description :  Convert gas proper mass density and pressure to sound speed squared
//
// Note        :  Eq. 14 in "Tseng et al. 2021, MNRAS, 504, 3298"
//
// Parameter   :  Dens       : Gas proper mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Sound speed squared
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_TaubMathews( const real Dens, const real Pres, const real Passive[],
                                           const double AuxArray_Flt[], const int AuxArray_Int[],
                                           const real *const Table[EOS_NTABLE_MAX] )
{

   real Cs2, Temp, factor;

   Temp   = Pres/Dens;
   factor = SQRT( (real)2.25*Temp*Temp + (real)1.0 );
   Cs2    = (real) 4.5*SQR(Temp) + (real) 5.0*Temp*factor;
   Cs2   /= (real)18.0*SQR(Temp) + (real)12.0*Temp*factor + (real)3.0;

#  ifdef GAMER_DEBUG
   if ( Cs2 >= (real)1.0  ||  Cs2 < (real)0.0 )
      printf( "ERROR : incorrect sound speed squared %14.7e (Dens %13.7e, Pres %13.7e) in %s !!\n",
              Cs2, Dens, Pres, __FUNCTION__ );
#  endif

   return Cs2;

} // FUNCTION : EoS_DensPres2CSqr_TaubMathews



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_GUESS_t  EoS_GuessHTilde_Ptr   = EoS_GuessHTilde_TaubMathews;
FUNC_SPACE EoS_H2TEM_t  EoS_HTilde2Temp_Ptr   = EoS_HTilde2Temp_TaubMathews;
FUNC_SPACE EoS_TEM2H_t  EoS_Temp2HTilde_Ptr   = EoS_Temp2HTilde_TaubMathews;
FUNC_SPACE EoS_DP2C_t   EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_TaubMathews;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_SetCPU/GPUFunc_TaubMathews
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_TaubMathews()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_SetFunc_TaubMathews( CPU_FuncPtr, GPU_FuncPtr );
//
//                3. Call-by-reference
//
// Parameter   :  EoS_HTilde2Temp_CPU/GPUPtr   : CPU/GPU function pointers to be set
//                EoS_Temp2HTilde_CPU/GPUPtr   : ...
//                EoS_DensPres2CSqr_CPU/GPUPtr : ...
//
// Return      :  EoS_HTilde2Temp_CPU, EoS_Temp2HTilde_CPU/GPUPtr, EoS_DensPres2CSqr_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void EoS_SetGPUFunc_TaubMathews( EoS_GUESS_t &EoS_GuessHTilde_GPUPtr,
                                 EoS_H2TEM_t &EoS_HTilde2Temp_GPUPtr,
                                 EoS_TEM2H_t &EoS_Temp2HTilde_GPUPtr,
                                 EoS_DP2C_t  &EoS_DensPres2CSqr_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_GuessHTilde_GPUPtr,   EoS_GuessHTilde_Ptr,   sizeof(EoS_GUESS_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_HTilde2Temp_GPUPtr,   EoS_HTilde2Temp_Ptr,   sizeof(EoS_H2TEM_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_Temp2HTilde_GPUPtr,   EoS_Temp2HTilde_Ptr,   sizeof(EoS_TEM2H_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t ) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_TaubMathews( EoS_GUESS_t &EoS_GuessHTilde_CPUPtr,
                                 EoS_H2TEM_t &EoS_HTilde2Temp_CPUPtr,
                                 EoS_TEM2H_t &EoS_Temp2HTilde_CPUPtr,
                                 EoS_DP2C_t  &EoS_DensPres2CSqr_CPUPtr )
{
   EoS_GuessHTilde_CPUPtr   = EoS_GuessHTilde_Ptr;
   EoS_HTilde2Temp_CPUPtr   = EoS_HTilde2Temp_Ptr;
   EoS_Temp2HTilde_CPUPtr   = EoS_Temp2HTilde_Ptr;
   EoS_DensPres2CSqr_CPUPtr = EoS_DensPres2CSqr_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void EoS_SetAuxArray_TaubMathews( double [] );
void EoS_SetCPUFunc_TaubMathews(EoS_GUESS_t &, EoS_H2TEM_t &, EoS_TEM2H_t &, EoS_DP2C_t & );
#ifdef GPU
void EoS_SetGPUFunc_TaubMathews(EoS_GUESS_t &, EoS_H2TEM_t &, EoS_TEM2H_t &, EoS_DP2C_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_TaubMathews
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
void EoS_Init_TaubMathews()
{

   EoS_SetAuxArray_TaubMathews( EoS_AuxArray_Flt, EoS_AuxArray_Int );
   EoS_SetCPUFunc_TaubMathews( EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_Temp2HTilde_CPUPtr, EoS_DensPres2CSqr_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_TaubMathews( EoS_GuessHTilde_GPUPtr, EoS_HTilde2Temp_GPUPtr, EoS_Temp2HTilde_GPUPtr, EoS_DensPres2CSqr_GPUPtr );
#  endif

} // FUNCTION : EoS_Init_TaubMathews

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO && defined SRHD )
