#include "Global.h"

#ifdef EXACT_COOLING


// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "CUFLU.h"
#include "CUDA_CheckError.h"
#include "CUFLU_Shared_FluUtility.cu"
#include "CUDA_ConstMemory.h"
#ifdef DUAL_ENERGY
#include "CUFLU_Shared_DualEnergy.cu"
#endif

extern double *d_SrcEC_TEF_lambda;
extern double *d_SrcEC_TEF_alpha;
extern double *d_SrcEC_TEFc;

#endif // #ifdef __CUDACC__


// local function prototypes
#ifndef __CUDACC__

void Src_SetAuxArray_ExactCooling( double [], int [] );
void Src_SetConstMemory_ExactCooling( const double AuxArray_Flt[], const int AuxArray_Int[],
                                      double *&DevPtr_Flt, int *&DevPtr_Int );
void Src_PassData2GPU_ExactCooling();
void Src_SetCPUFunc_ExactCooling( SrcFunc_t & );
#ifdef GPU
void Src_SetGPUFunc_ExactCooling( SrcFunc_t & );
#endif
#ifdef GPU
void CUAPI_MemFree_ExactCooling();
#endif
void Src_WorkBeforeMajorFunc_ExactCooling( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                           double AuxArray_Flt[], int AuxArray_Int[] );
void Src_End_ExactCooling();

void Cool_fct( double Dens, double Temp, double* Emis, double* Lambdat, double Z, double cl_moli_mole, double mp );
#endif // #ifndef __CUDACC__
GPU_DEVICE static
double TEF( double TEMP, int k, const double TEF_lambda[], const double TEF_alpha[], const double TEFc[],
            const double AuxArray_Flt[], const int AuxArray_Int[] );
GPU_DEVICE static
double TEFinv( double Y, int k, const double TEF_lambda[], const double TEF_alpha[], const double TEFc[],
               const double AuxArray_Flt[], const int AuxArray_Int[] );


/********************************************************
1. Exact-cooling source term
   --> Enabled by the runtime option "SRC_EXACTCOOLING"

2. This file is shared by both CPU and GPU

   CUSRC_Src_ExactCooling.cu -> CPU_Src_ExactCooling.cpp

3. Four steps are required to implement a source term

   I.   Set auxiliary arrays
   II.  Implement the source-term function
   III. [Optional] Add the work to be done every time
        before calling the major source-term function
   IV.  Set initialization functions

4. The source-term function must be thread-safe and
   not use any global variable
********************************************************/



// =======================
// I. Set auxiliary arrays
// =======================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetAuxArray_ExactCooling
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by Src_Init_ExactCooling()
//                2. AuxArray_Flt/Int[] have the size of SRC_NAUX_EC defined in Macro.h (default = 10)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_SetAuxArray_ExactCooling( double AuxArray_Flt[], int AuxArray_Int[] )
{

   const int    TEF_N      = SrcTerms.EC_TEF_N;      // number of points for lambda(T) sampling in LOG
   const bool   subcycling = SrcTerms.EC_subcycling; // whether to use subcycling
   const int    TEF_int    = TEF_N-1;                // number of intervals
   const double TEF_TN     = 1.e14;                  // == Tref, must be high enough, but affects sampling resolution (Kelvin)
   const double TEF_Tmin   = MIN_TEMP;               // MIN temperature
   if ( TEF_Tmin <= 0.0 )
      Aux_Error( ERROR_INFO, "TEF_Tmin (%14.7e) <= 0.0 !!\n", TEF_Tmin );
   const double TEF_dltemp   = (log10(TEF_TN) - log10(TEF_Tmin))/TEF_int; // sampling resolution (Kelvin), LOG!

   const double cl_X         = 0.7;                                       // mass-fraction of hydrogen
   const double cl_Z         = 0.018;                                     // metallicity (in Zsun)
   const double cl_mol       = 1.0/(2*cl_X+0.75*(1-cl_X-cl_Z)+cl_Z*0.5);  // mean (total)  molecular weights
   const double cl_mole      = 2.0/(1+cl_X);                              // mean electron molecular weights
   const double cl_moli      = 1.0/cl_X;                                  // mean proton   molecular weights
   const double cl_moli_mole = cl_moli*cl_mole;

// store them in the aux array
   AuxArray_Flt[0] = 1.0/(GAMMA-1.0);
   AuxArray_Flt[1] = TEF_TN;
   AuxArray_Flt[2] = TEF_Tmin;
   AuxArray_Flt[3] = TEF_dltemp;
   AuxArray_Flt[4] = cl_Z;
   AuxArray_Flt[5] = cl_moli_mole;
   AuxArray_Flt[6] = cl_mol;
   AuxArray_Flt[7] = MU_NORM/UNIT_M;                       // mp: proton mass
   AuxArray_Flt[8] = (Const_kB/UNIT_E) * (MU_NORM/UNIT_M); // kB*mp

   AuxArray_Int[0] = TEF_N;
   AuxArray_Int[1] = subcycling;

} // FUNCTION : Src_SetAuxArray_ExactCooling
#endif // #ifndef __CUDACC__


// ======================================
// II. Implement the source-term function
// ======================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_ExactCooling
// Description :  Major source-term function
//
// Note        :  1. Invoked by CPU/GPU_SrcSolver_IterateAllCells()
//                2. See Src_SetAuxArray_ExactCooling() for the values stored in AuxArray_Flt/Int[]
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
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                EoS               : EoS object
//                AuxArray_*        : Auxiliary arrays (see the Note above)
//
// Return      :  fluid[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void Src_ExactCooling( real fluid[], const real B[], const SrcTerms_t *SrcTerms, const real dt,
                              const real dh, const double x, const double y, const double z,
                              const double TimeNew, const double TimeOld, const real MinDens,
                              const real MinPres, const real MinEint, const long PassiveFloor, const EoS_t *EoS,
                              const double AuxArray_Flt[], const int AuxArray_Int[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );
#  endif

} // FUNCTION : Src_ExactCooling



//-------------------------------------------------------------------------------------------------------
// Function    :  TEF
// Description :  the temporal evolution function (TEF)
//
// Note        :  TBF
//
// Parameter   :  TBF
//
// Return      :  TBF
//-----------------------------------------------------------------------------------------
GPU_DEVICE static
double TEF( double TEMP, int k, const double TEF_lambda[], const double TEF_alpha[], const double TEFc[],
            const double AuxArray_Flt[], const int AuxArray_Int[] )
{
   return 0.0;
} // FUNCTION : TEF



//-------------------------------------------------------------------------------------------------------
// Function    :  TEFinv
// Description :  the INVERSE temporal evolution function (TEF^-1)
//
// Note        :  TBF
//
// Parameter   :  TBF
//
// Return      :  TBF
//-----------------------------------------------------------------------------------------
GPU_DEVICE static
double TEFinv( double Y, int k, const double TEF_lambda[], const double TEF_alpha[], const double TEFc[],
               const double AuxArray_Flt[], const int AuxArray_Int[] )
{
   return 0.0;
} // FUNCTION : TEFinv



// ==================================================
// III. [Optional] Add the work to be done every time
//      before calling the major source-term function
// ==================================================

#ifndef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_WorkBeforeMajorFunc_ExactCooling
// Description :  Specify work to be done every time before calling the major source-term function
//
// Note        :  1. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  lv               : Target refinement level
//                TimeNew          : Target physical time to reach
//                TimeOld          : Physical time before update
//                                   --> The major source-term function will update the system from TimeOld to TimeNew
//                dt               : Time interval to advance solution
//                                   --> Physical coordinates : TimeNew - TimeOld == dt
//                                       Comoving coordinates : TimeNew - TimeOld == delta(scale factor) != dt
//                AuxArray_Flt/Int : Auxiliary arrays
//                                   --> Can be used and/or modified here
//                                   --> Must call Src_SetConstMemory_ExactCooling() after modification
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void Src_WorkBeforeMajorFunc_ExactCooling( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                           double AuxArray_Flt[], int AuxArray_Int[] )
{
//  nothing to do here
} // FUNCTION : Src_WorkBeforeMajorFunc_ExactCooling
#endif // #ifndef __CUDACC__



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_PassData2GPU_ExactCooling
// Description :  Transfer data to GPU
//
// Note        :  1. Invoked by Src_Init_ExactCooling()
//                2. Use synchronous transfer
//
// Parameter   :  None
// Return      :  None
// -------------------------------------------------------------------------------------------------------
void Src_PassData2GPU_ExactCooling()
{
   const long EC_TEF_MemSize = sizeof(double)*SrcTerms.EC_TEF_N;

   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_SrcEC_TEF_lambda, EC_TEF_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_SrcEC_TEF_alpha,  EC_TEF_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_SrcEC_TEFc,       EC_TEF_MemSize )  );

// store the device pointers in SrcTerms when using GPU
   SrcTerms.EC_TEF_lambda_DevPtr = d_SrcEC_TEF_lambda;
   SrcTerms.EC_TEF_alpha_DevPtr  = d_SrcEC_TEF_alpha;
   SrcTerms.EC_TEFc_DevPtr       = d_SrcEC_TEFc;

// use synchronous transfer
   CUDA_CHECK_ERROR(  cudaMemcpy( d_SrcEC_TEF_lambda, h_SrcEC_TEF_lambda, EC_TEF_MemSize, cudaMemcpyHostToDevice )  );
   CUDA_CHECK_ERROR(  cudaMemcpy( d_SrcEC_TEF_alpha,  h_SrcEC_TEF_alpha,  EC_TEF_MemSize, cudaMemcpyHostToDevice )  );
   CUDA_CHECK_ERROR(  cudaMemcpy( d_SrcEC_TEFc,       h_SrcEC_TEFc,       EC_TEF_MemSize, cudaMemcpyHostToDevice )  );

} // FUNCTION : Src_PassData2GPU_ExactCooling
#endif // #ifdef __CUDACC__


// ================================
// IV. Set initialization functions
// ================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE SrcFunc_t SrcFunc_Ptr = Src_ExactCooling;

//-----------------------------------------------------------------------------------------
// Function    :  Src_SetCPU/GPUFunc_ExactCooling
// Description :  Return the function pointer of the CPU/GPU source-term function
//
// Note        :  1. Invoked by Src_Init_ExactCooling()
//                2. Call-by-reference
//
// Parameter   :  SrcFunc_CPU/GPUPtr : CPU/GPU function pointer to be set
//
// Return      :  SrcFunc_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void Src_SetGPUFunc_ExactCooling( SrcFunc_t &SrcFunc_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &SrcFunc_GPUPtr, SrcFunc_Ptr, sizeof(SrcFunc_t) )  );
}

#else

void Src_SetCPUFunc_ExactCooling( SrcFunc_t &SrcFunc_CPUPtr )
{
   SrcFunc_CPUPtr = SrcFunc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetConstMemory_ExactCooling
// Description :  Set the constant memory variables on GPU
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by Src_Init_ExactCooling() and, if necessary, Src_WorkBeforeMajorFunc_ExactCooling()
//                3. SRC_NAUX_EC is defined in Macro.h
//
// Parameter   :  AuxArray_Flt/Int : Auxiliary arrays to be copied to the constant memory
//                DevPtr_Flt/Int   : Pointers to store the addresses of constant memory arrays
//
// Return      :  c_Src_EC_AuxArray_Flt[], c_Src_EC_AuxArray_Int[], DevPtr_Flt, DevPtr_Int
//---------------------------------------------------------------------------------------------------
void Src_SetConstMemory_ExactCooling( const double AuxArray_Flt[], const int AuxArray_Int[],
                                      double *&DevPtr_Flt, int *&DevPtr_Int )
{

// copy data to constant memory
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_EC_AuxArray_Flt, AuxArray_Flt, SRC_NAUX_EC*sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_EC_AuxArray_Int, AuxArray_Int, SRC_NAUX_EC*sizeof(int   ) )  );

// obtain the constant-memory pointers
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Flt, c_Src_EC_AuxArray_Flt )  );
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Int, c_Src_EC_AuxArray_Int )  );

} // FUNCTION : Src_SetConstMemory_ExactCooling
#endif // #ifdef __CUDACC__



#ifndef __CUDACC__
//-----------------------------------------------------------------------------------------
// Function    :  Src_Init_ExactCooling
// Description :  Initialize the exact-cooling source term
//
// Note        :  1. Set auxiliary arrays by invoking Src_SetAuxArray_*()
//                   --> Copy to the GPU constant memory and store the associated addresses
//                2. Set the source-term function by invoking Src_SetCPU/GPUFunc_*()
//                3. Invoked by Src_Init()
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_Init_ExactCooling()
{

   Aux_Error( ERROR_INFO, "SRC_EXACTCOOLING is not supported !!\n" );

// set the auxiliary arrays
   Src_SetAuxArray_ExactCooling( Src_EC_AuxArray_Flt, Src_EC_AuxArray_Int );

// copy the auxiliary arrays to the GPU constant memory and store the associated addresses
#  ifdef GPU
   Src_SetConstMemory_ExactCooling( Src_EC_AuxArray_Flt, Src_EC_AuxArray_Int,
                                    SrcTerms.EC_AuxArrayDevPtr_Flt, SrcTerms.EC_AuxArrayDevPtr_Int );
#  else
   SrcTerms.EC_AuxArrayDevPtr_Flt = Src_EC_AuxArray_Flt;
   SrcTerms.EC_AuxArrayDevPtr_Int = Src_EC_AuxArray_Int;
#  endif

// set the major source-term function
   Src_SetCPUFunc_ExactCooling( SrcTerms.EC_CPUPtr );

#  ifdef GPU
   Src_SetGPUFunc_ExactCooling( SrcTerms.EC_GPUPtr );
   SrcTerms.EC_FuncPtr = SrcTerms.EC_GPUPtr;
#  else
   SrcTerms.EC_FuncPtr = SrcTerms.EC_CPUPtr;
#  endif

   if ( OPT__INIT == INIT_BY_RESTART )
      for (int i=0; i<NLEVEL; i++)   IsInit_tcool[i] = true;
   else
      for (int i=0; i<NLEVEL; i++)   IsInit_tcool[i] = false;

// allocate h_SrcEC_* arrays
   h_SrcEC_TEF_lambda = new double [SrcTerms.EC_TEF_N];
   h_SrcEC_TEF_alpha  = new double [SrcTerms.EC_TEF_N];
   h_SrcEC_TEFc       = new double [SrcTerms.EC_TEF_N];

#  ifdef GPU
   Src_PassData2GPU_ExactCooling();
#  else
   SrcTerms.EC_TEF_lambda_DevPtr = h_SrcEC_TEF_lambda;
   SrcTerms.EC_TEF_alpha_DevPtr  = h_SrcEC_TEF_alpha;
   SrcTerms.EC_TEFc_DevPtr       = h_SrcEC_TEFc;
#  endif

} // FUNCTION : Src_Init_ExactCooling



//-------------------------------------------------------------------------------------------------------
// Function    :  Cool_fct
// Description :  Sutherland-Dopita cooling function, with optimal parmetrization over a wide range of T and Z
//
// Note        :  TBF
//
// Parameter   :  TBF
//
// Return      :  TBF
//-----------------------------------------------------------------------------------------
void Cool_fct( double Dens, double Temp, double* Emis, double* Lambdat, double Z, double cl_moli_mole, double mp )
{
} // FUNCTION : Cool_fct



//-----------------------------------------------------------------------------------------
// Function    :  Src_End_ExactCooling
// Description :  Free the resources used by the exact-cooling source term
//
// Note        :  1. Invoked by Src_End()
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_End_ExactCooling()
{

// free the memory of the h_SrcEC_* arrays
   delete [] h_SrcEC_TEF_lambda;     h_SrcEC_TEF_lambda = NULL;
   delete [] h_SrcEC_TEF_alpha;      h_SrcEC_TEF_alpha  = NULL;
   delete [] h_SrcEC_TEFc;           h_SrcEC_TEFc       = NULL;

#  ifdef GPU
   CUAPI_MemFree_ExactCooling();
#  else
   SrcTerms.EC_TEF_lambda_DevPtr = NULL;
   SrcTerms.EC_TEF_alpha_DevPtr  = NULL;
   SrcTerms.EC_TEFc_DevPtr       = NULL;
#  endif

} // FUNCTION : Src_End_ExactCooling
#endif // #ifndef __CUDACC__



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemFree_ExactCooling
// Description :  Free the GPU memory of the ExactCooling arrays
//
// Note        :  1. Invoked by Src_End_ExactCooling()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemFree_ExactCooling()
{

   if ( d_SrcEC_TEF_lambda != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_SrcEC_TEF_lambda )  );  d_SrcEC_TEF_lambda = NULL; }
   if ( d_SrcEC_TEF_alpha  != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_SrcEC_TEF_alpha  )  );  d_SrcEC_TEF_alpha  = NULL; }
   if ( d_SrcEC_TEFc       != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_SrcEC_TEFc       )  );  d_SrcEC_TEFc       = NULL; }

   SrcTerms.EC_TEF_lambda_DevPtr = NULL;
   SrcTerms.EC_TEF_alpha_DevPtr  = NULL;
   SrcTerms.EC_TEFc_DevPtr       = NULL;

} // FUNCTION : CUAPI_MemFree_ExactCooling
#endif // #ifdef __CUDACC__


#endif // #ifdef EXACT_COOLING
