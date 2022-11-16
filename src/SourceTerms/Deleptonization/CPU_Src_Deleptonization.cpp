#include "CUFLU.h"

#if ( MODEL == HYDRO )



// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "Global.h"
#include "CUDA_CheckError.h"
#include "CUFLU_Shared_FluUtility.cu"
#include "CUDA_ConstMemory.h"

extern real (*d_SrcDlepProf_Data)[SRC_DLEP_PROF_NBINMAX];
extern real  *d_SrcDlepProf_Radius;

#endif // #ifdef __CUDACC__


// local function prototypes
#ifndef __CUDACC__

void Src_SetAuxArray_Deleptonization( double [], int [] );
void Src_SetCPUFunc_Deleptonization( SrcFunc_t & );
#ifdef GPU
void Src_SetGPUFunc_Deleptonization( SrcFunc_t & );
#endif
void Src_SetConstMemory_Deleptonization( const double AuxArray_Flt[], const int AuxArray_Int[],
                                         double *&DevPtr_Flt, int *&DevPtr_Int );
void Src_PassData2GPU_Deleptonization();

#endif



/********************************************************
1. Deleptonization source term
   --> Enabled by the runtime option "SRC_DELEPTONIZATION"

2. This file is shared by both CPU and GPU

   CUSRC_Src_Deleptonization.cu -> CPU_Src_Deleptonization.cpp

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
// Function    :  Src_SetAuxArray_Deleptonization
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by Src_Init_Deleptonization()
//                2. AuxArray_Flt/Int[] have the size of SRC_NAUX_DLEP defined in Macro.h (default = 5)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_SetAuxArray_Deleptonization( double AuxArray_Flt[], int AuxArray_Int[] )
{

// TBF

} // FUNCTION : Src_SetAuxArray_Deleptonization
#endif // #ifndef __CUDACC__



// ======================================
// II. Implement the source-term function
// ======================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Deleptonization
// Description :  Major source-term function
//
// Note        :  1. Invoked by CPU/GPU_SrcSolver_IterateAllCells()
//                2. See Src_SetAuxArray_Deleptonization() for the values stored in AuxArray_Flt/Int[]
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
//                EoS               : EoS object
//                AuxArray_*        : Auxiliary arrays (see the Note above)
//
// Return      :  fluid[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void Src_Deleptonization( real fluid[], const real B[],
                                 const SrcTerms_t *SrcTerms, const real dt, const real dh,
                                 const double x, const double y, const double z,
                                 const double TimeNew, const double TimeOld,
                                 const real MinDens, const real MinPres, const real MinEint,
                                 const EoS_t *EoS, const double AuxArray_Flt[], const int AuxArray_Int[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );
#  endif

// TBF
// profiles are stored in SrcTerms->Dlep_Profile_DataDevPtr/Dlep_Profile_RadiusDevPtr/Dlep_Profile_NBin
// --> see "include/SrcTerms.h"

} // FUNCTION : Src_Deleptonization



// ==================================================
// III. [Optional] Add the work to be done every time
//      before calling the major source-term function
// ==================================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_WorkBeforeMajorFunc_Deleptonization
// Description :  Specify work to be done every time before calling the major source-term function
//
// Note        :  1. Invoked by Src_WorkBeforeMajorFunc()
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
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
//                                   --> Must call Src_SetConstMemory_Deleptonization() after modification
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_WorkBeforeMajorFunc_Deleptonization( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                              double AuxArray_Flt[], int AuxArray_Int[] )
{

// TBF

/*
// compute profiles
// --> here is just an example; see GREP for a more efficient implementation
// --> SRC_DLEP_PROF_NVAR and SRC_DLEP_PROF_NBINMAX are defined in Macro.h (default = 6 and 4000, respectively)
// --> be careful about the issue of center drifting
   const double      Center[3]      = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double      MaxRadius      = 0.5*amr->BoxSize[0];
   const double      MinBinSize     = amr->dh[MAX_LEVEL];
   const bool        LogBin         = true;
   const double      LogBinRatio    = 1.25;
   const bool        RemoveEmptyBin = true;
   const long        TVar[]         = { _DENS, _MOMX, _ENGY, _PRES, _VELR, _EINT };
   const int         NProf          = SRC_DLEP_PROF_NVAR;
   const int         SingleLv       = -1;
   const int         MaxLv          = -1;
   const PatchType_t PatchType      = PATCH_LEAF;
   const double      PrepTime       = TimeNew;

   Profile_t *Prof[SRC_DLEP_PROF_NVAR];
   for (int v=0; v<SRC_DLEP_PROF_NVAR; v++)  Prof[v] = new Profile_t();

   Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin,
                       TVar, NProf, SingleLv, MaxLv, PatchType, PrepTime );


// check and store the number of radial bins
   if ( Prof[0]->NBin > SRC_DLEP_PROF_NBINMAX )
      Aux_Error( ERROR_INFO, "Number of radial bins (%d) exceeds the maximum size (%d) !!\n",
                 Prof[0]->NBin, SRC_DLEP_PROF_NBINMAX );

   SrcTerms.Dlep_Profile_NBin = Prof[0]->NBin;


// store profiles in the host arrays
// --> note the typecasting from double to real
   for (int v=0; v<SRC_DLEP_PROF_NVAR; v++)
   for (int b=0; b<Prof[v]->NBin; b++)
      h_SrcDlepProf_Data[v][b] = (real)Prof[v]->Data[b];

   for (int b=0; b<Prof[0]->NBin; b++)
      h_SrcDlepProf_Radius[b] = (real)Prof[0]->Radius[b];


// pass profiles to GPU
#  ifdef GPU
   Src_PassData2GPU_Deleptonization();
#  endif


// uncomment the following lines if the auxiliary arrays have been modified
//#  ifdef GPU
//   Src_SetConstMemory_Deleptonization( AuxArray_Flt, AuxArray_Int,
//                                       SrcTerms.Dlep_AuxArrayDevPtr_Flt, SrcTerms.Dlep_AuxArrayDevPtr_Int );
//#  endif


// free memory
   for (int v=0; v<SRC_DLEP_PROF_NVAR; v++)  delete Prof[v];
*/

} // FUNCTION : Src_WorkBeforeMajorFunc_Deleptonization
#endif



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_PassData2GPU_Deleptonization
// Description :  Transfer data to GPU
//
// Note        :  1. Invoked by Src_WorkBeforeMajorFunc_Deleptonization()
//                2. Use synchronous transfer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Src_PassData2GPU_Deleptonization()
{

   const long Size_Data   = sizeof(real)*SRC_DLEP_PROF_NVAR*SRC_DLEP_PROF_NBINMAX;
   const long Size_Radius = sizeof(real)*                   SRC_DLEP_PROF_NBINMAX;

// use synchronous transfer
   CUDA_CHECK_ERROR(  cudaMemcpy( d_SrcDlepProf_Data,   h_SrcDlepProf_Data,   Size_Data,   cudaMemcpyHostToDevice )  );
   CUDA_CHECK_ERROR(  cudaMemcpy( d_SrcDlepProf_Radius, h_SrcDlepProf_Radius, Size_Radius, cudaMemcpyHostToDevice )  );

} // FUNCTION : Src_PassData2GPU_Deleptonization
#endif // #ifdef __CUDACC__



// ================================
// IV. Set initialization functions
// ================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE SrcFunc_t SrcFunc_Ptr = Src_Deleptonization;

//-----------------------------------------------------------------------------------------
// Function    :  Src_SetCPU/GPUFunc_Deleptonization
// Description :  Return the function pointer of the CPU/GPU source-term function
//
// Note        :  1. Invoked by Src_Init_Deleptonization()
//                2. Call-by-reference
//
// Parameter   :  SrcFunc_CPU/GPUPtr : CPU/GPU function pointer to be set
//
// Return      :  SrcFunc_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void Src_SetGPUFunc_Deleptonization( SrcFunc_t &SrcFunc_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &SrcFunc_GPUPtr, SrcFunc_Ptr, sizeof(SrcFunc_t) )  );
}

#else

void Src_SetCPUFunc_Deleptonization( SrcFunc_t &SrcFunc_CPUPtr )
{
   SrcFunc_CPUPtr = SrcFunc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetConstMemory_Deleptonization
// Description :  Set the constant memory variables on GPU
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by Src_Init_Deleptonizatio() and, if necessary, Src_WorkBeforeMajorFunc_Deleptonizatio()
//                3. SRC_NAUX_DLEP is defined in Macro.h
//
// Parameter   :  AuxArray_Flt/Int : Auxiliary arrays to be copied to the constant memory
//                DevPtr_Flt/Int   : Pointers to store the addresses of constant memory arrays
//
// Return      :  c_Src_Dlep_AuxArray_Flt[], c_Src_Dlep_AuxArray_Int[], DevPtr_Flt, DevPtr_Int
//---------------------------------------------------------------------------------------------------
void Src_SetConstMemory_Deleptonization( const double AuxArray_Flt[], const int AuxArray_Int[],
                                         double *&DevPtr_Flt, int *&DevPtr_Int )
{

// copy data to constant memory
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_Dlep_AuxArray_Flt, AuxArray_Flt, SRC_NAUX_DLEP*sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_Dlep_AuxArray_Int, AuxArray_Int, SRC_NAUX_DLEP*sizeof(int   ) )  );

// obtain the constant-memory pointers
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Flt, c_Src_Dlep_AuxArray_Flt )  );
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Int, c_Src_Dlep_AuxArray_Int )  );

} // FUNCTION : Src_SetConstMemory_Deleptonization
#endif // #ifdef __CUDACC__



#ifndef __CUDACC__

//-----------------------------------------------------------------------------------------
// Function    :  Src_Init_Deleptonization
// Description :  Initialize the deleptonization source term
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
void Src_Init_Deleptonization()
{

// set the auxiliary arrays
   Src_SetAuxArray_Deleptonization( Src_Dlep_AuxArray_Flt, Src_Dlep_AuxArray_Int );

// copy the auxiliary arrays to the GPU constant memory and store the associated addresses
#  ifdef GPU
   Src_SetConstMemory_Deleptonization( Src_Dlep_AuxArray_Flt, Src_Dlep_AuxArray_Int,
                                       SrcTerms.Dlep_AuxArrayDevPtr_Flt, SrcTerms.Dlep_AuxArrayDevPtr_Int );
#  else
   SrcTerms.Dlep_AuxArrayDevPtr_Flt = Src_Dlep_AuxArray_Flt;
   SrcTerms.Dlep_AuxArrayDevPtr_Int = Src_Dlep_AuxArray_Int;
#  endif

// set the major source-term function
   Src_SetCPUFunc_Deleptonization( SrcTerms.Dlep_CPUPtr );

#  ifdef GPU
   Src_SetGPUFunc_Deleptonization( SrcTerms.Dlep_GPUPtr );
   SrcTerms.Dlep_FuncPtr = SrcTerms.Dlep_GPUPtr;
#  else
   SrcTerms.Dlep_FuncPtr = SrcTerms.Dlep_CPUPtr;
#  endif

} // FUNCTION : Src_Init_Deleptonization



//-----------------------------------------------------------------------------------------
// Function    :  Src_End_Deleptonization
// Description :  Release the resources used by the deleptonization source term
//
// Note        :  1. Invoked by Src_End()
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_End_Deleptonization()
{

// TBF

} // FUNCTION : Src_End_Deleptonization

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO )
