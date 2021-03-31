#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY




// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__

extern double IsothermalSlab_Center[3];
extern double IsothermalSlab_Temperature;
extern double IsothermalSlab_PeakDens;
extern double IsothermalSlab_Truncation;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_IsothermalSlab
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_IsothermalSlab()
//
// Note        :  1. Invoked by Init_ExtPot_IsothermalSlab()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_IsothermalSlab( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[ 0] = IsothermalSlab_Center[0];                // x coordinate of the IsothermalSlab center
   AuxArray_Flt[ 1] = IsothermalSlab_Center[1];                // y ...
   AuxArray_Flt[ 2] = IsothermalSlab_Center[2];                // z ...
   AuxArray_Flt[ 3] = IsothermalSlab_Temperature;        // 
   AuxArray_Flt[ 4] = IsothermalSlab_PeakDens;           // 
   AuxArray_Flt[ 6] = NEWTON_G;
   AuxArray_Flt[ 7] =  amr->BoxSize[2];
   AuxArray_Flt[ 8] = IsothermalSlab_Truncation;

} // FUNCTION : SetExtPotAuxArray_IsothermalSlab
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_IsothermalSlab
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_IsothermalSlab()
//
// Parameter   :  x/y/z             : Target spatial coordinates
//                Time              : Target physical time
//                UserArray_Flt/Int : User-provided floating-point/integer auxiliary arrays
//                Usage             : Different usages of external potential when computing total potential on level Lv
//                                    --> EXT_POT_USAGE_ADD     : add external potential on Lv
//                                        EXT_POT_USAGE_SUB     : subtract external potential for preparing self-gravity potential on Lv-1
//                                        EXT_POT_USAGE_SUB_TINT: like SUB but for temporal interpolation
//                                    --> This parameter is useless in most cases
//                PotTable          : 3D potential table used by EXT_POT_TABLE
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_IsothermalSlab( const double x, const double y, const double z, const double Time,
                             const double UserArray_Flt[], const int UserArray_Int[],
                             const ExtPotUsage_t Usage, const real PotTable[] )
{

// halo potential
   const double cz                                =       UserArray_Flt[ 2];   // z ...
   const real   IsothermalSlab_Temperature        = (real)UserArray_Flt[ 3];   // 
   const real   IsothermalSlab_PeakDens           = (real)UserArray_Flt[ 4];   // 
   const real   NewtonG                           =       UserArray_Flt[ 6];
   const real   BoxSize_Z                         =       UserArray_Flt[ 7];    
   const double dz                                = z - cz;
   const real   IsothermalSlab_Truncation         = (real)UserArray_Flt[ 8];   // 


   real stellarDiskPot = sqrt( ( 2.0*M_PI*NewtonG*IsothermalSlab_PeakDens ) / IsothermalSlab_Temperature );

#  ifndef __CUDACC__
   if ( cz != 0.5*BoxSize_Z ) Aux_Error( ERROR_INFO, "We expect the z-position of stellar disk is at box-center!\n"); 
#  endif

   if ( abs(dz) > IsothermalSlab_Truncation )
     stellarDiskPot = 2.0*IsothermalSlab_Temperature*log(cosh(IsothermalSlab_Truncation*stellarDiskPot));
   else
     stellarDiskPot = 2.0*IsothermalSlab_Temperature*log(cosh(dz*stellarDiskPot));


   return stellarDiskPot;

} // FUNCTION : ExtPot_IsothermalSlab



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_IsothermalSlab;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_IsothermalSlab
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_IsothermalSlab()
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_IsothermalSlab( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_IsothermalSlab( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_IsothermalSlab( double [], int [] );
void SetCPUExtPot_IsothermalSlab( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_IsothermalSlab( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_IsothermalSlab
// Description :  Initialize external potential
//
// Note        :  1. Set auxiliary arrays by invoking SetExtPotAuxArray_*()
//                   --> They will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external potential major routines by invoking SetCPU/GPUExtPot_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtPot_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtPot_IsothermalSlab()
{

   SetExtPotAuxArray_IsothermalSlab( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int );
   SetCPUExtPot_IsothermalSlab( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_IsothermalSlab( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_IsothermalSlab

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
