#include "CUFLU.h"



// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "CUAPI.h"
#include "CUFLU_Shared_FluUtility.cu"
#include "CUDA_ConstMemory.h"

#endif // #ifdef __CUDACC__



/********************************************************
1. Template of a user-defined source term
   --> Enabled by the runtime option "SRC_USER"

2. This file is shared by both CPU and GPU

   CUSRC_Src_User_Template.cu -> CPU_Src_User_Template.cpp

3. Three steps are required to implement a source term

   I.   Set auxiliary arrays
   II.  Implement the source-term function
   III. Set initialization functions

4. The source-term function must be thread-safe and
   not use any global variable
********************************************************/



// =============================================
// I. Set auxiliary arrays
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetAuxArray_User_Template
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by Src_Init_User_Template()
//                2. AuxArray_Flt/Int[] have the size of SRC_NAUX_USER defined in Macro.h (default = 10)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_SetAuxArray_User_Template( double AuxArray_Flt[], int AuxArray_Int[] )
{

   /*
   AuxArray_Flt[0] = ...;
   AuxArray_Flt[1] = ...;

   AuxArray_Int[0] = ...;
   AuxArray_Int[1] = ...;
   */

} // FUNCTION : Src_SetAuxArray_User_Template
#endif // #ifndef __CUDACC__



// =============================================
// II. Implement the source-term function
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_User_Template
// Description :  Major source-term function
//
// Note        :  1. Invoked by CPU/GPU_SrcSolver_IterateAllCells()
//                2. See Src_SetAuxArray_User_Template() for the values stored in AuxArray_Flt/Int[]
//                3. Shared by both CPU and GPU
//
// Parameter   :  fluid             : Fluid array storing both the input and updated values
//                                    --> Including both active and passive variables
//                B                 : Cell-centered magnetic field
//                SrcTerms          : Structure storing all source-term variables
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                x/y/z             : Target physical coordinates
//                TimeNew           : Target physical time to reach
//                TimeOld           : Physical time before update
//                                    --> This function updates physical time from TimeOld to TimeNew
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_DensPres2Eint : EoS routine to compute the gas internal energy
//                EoS_DensPres2CSqr : EoS routine to compute the sound speed square
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                EoS_Table         : EoS tables
//                AuxArray_*        : Auxiliary arrays (see the Note above)
//
// Return      :  fluid[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void Src_User_Template( real fluid[], const real B[],
                               const SrcTerms_t SrcTerms, const real dt, const real dh,
                               const double x, const double y, const double z,
                               const double TimeNew, const double TimeOld,
                               const real MinDens, const real MinPres, const real MinEint,
                               const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2E_t EoS_DensPres2Eint,
                               const EoS_DP2C_t EoS_DensPres2CSqr,
                               const double EoS_AuxArray_Flt[],
                               const int    EoS_AuxArray_Int[],
                               const real *const EoS_Table[EOS_NTABLE_MAX],
                               const double AuxArray_Flt[], const int AuxArray_Int[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );
#  endif

// example
   /*
   const bool CheckMinEint_Yes = true;
   const real CoolingRate      = (real)AuxArray_Flt[0];

   real Eint, Enth, Emag;

#  ifdef MHD
   Emag  = (real)0.5*( SQR(B[MAGX]) + SQR(B[MAGY]) + SQR(B[MAGZ]) );
#  else
   Emag  = (real)0.0;
#  endif
   Eint  = Hydro_Con2Eint( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                           CheckMinEint_Yes, MinEint, Emag );
   Enth  = fluid[ENGY] - Eint;
   Eint -= fluid[DENS]*CoolingRate*dt;

   fluid[ENGY] = Enth + Eint;
   */

} // FUNCTION : Src_User_Template



// =============================================
// III. Set initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE SrcFunc_t SrcFunc_Ptr = Src_User_Template;

//-----------------------------------------------------------------------------------------
// Function    :  Src_SetFunc_User_Template
// Description :  Return the function pointer of the CPU/GPU source-term function
//
// Note        :  1. Invoked by Src_Init_User_Template()
//                2. Call-by-reference
//                3. Use either CPU or GPU but not both of them
//
// Parameter   :  SrcFunc_CPU/GPUPtr : CPU/GPU function pointer to be set
//
// Return      :  SrcFunc_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void Src_SetFunc_User_Template( SrcFunc_t &SrcFunc_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &SrcFunc_GPUPtr, SrcFunc_Ptr, sizeof(SrcFunc_t) )  );
}

#elif ( !defined GPU )

void Src_SetFunc_User_Template( SrcFunc_t &SrcFunc_CPUPtr )
{
   SrcFunc_CPUPtr = SrcFunc_Ptr;
}

#endif // #ifdef __CUDACC__ ... elif ...



#ifndef __CUDACC__

// local function prototypes
void Src_SetAuxArray_User_Template( double [], int [] );
void Src_SetFunc_User_Template( SrcFunc_t & );

//-----------------------------------------------------------------------------------------
// Function    :  Src_Init_User_Template
// Description :  Initialize a user-specified source term
//
// Note        :  1. Set auxiliary arrays by invoking Src_SetAuxArray_*()
//                2. Set the source-term function by invoking Src_SetFunc_*()
//                   --> Unlike other modules (e.g., EoS), here we use either CPU or GPU but not
//                       both of them
//                3. Invoked by Src_Init()
//                   --> Enable it by linking to the function pointer "Src_Init_User_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_Init_User_Template()
{

   Src_SetAuxArray_User_Template( SRC_TERMS.User_AuxArray_Flt, SRC_TERMS.User_AuxArray_Int );
   Src_SetFunc_User_Template( SRC_TERMS.User_FuncPtr );

} // FUNCTION : Src_Init_User_Template



//-----------------------------------------------------------------------------------------
// Function    :  Src_End_User_Template
// Description :  Free the resources used by a user-specified source term
//
// Note        :  1. Invoked by Src_End()
//                   --> Enable it by linking to the function pointer "Src_End_User_Ptr"
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_End_User_Template()
{


} // FUNCTION : Src_End_User_Template

#endif // #ifndef __CUDACC__
