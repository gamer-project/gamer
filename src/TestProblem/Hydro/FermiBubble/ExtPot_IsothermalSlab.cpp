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
extern real IsothermalSlab_VelocityDispersion;
extern real IsothermalSlab_PeakDens;
extern real IsothermalSlab_Truncation;
extern real interfaceHeight;

extern  int Jet_Ambient;
extern real distance_h;
extern real v_halo;

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
   AuxArray_Flt[ 3] = IsothermalSlab_VelocityDispersion;       //
   AuxArray_Flt[ 4] = IsothermalSlab_PeakDens;                 //
   AuxArray_Flt[ 5] = NEWTON_G;
   AuxArray_Flt[ 6] = amr->BoxSize[2];
   AuxArray_Flt[ 7] = IsothermalSlab_Truncation;
   AuxArray_Flt[ 8] = distance_h;
   AuxArray_Flt[ 9] = v_halo;
   AuxArray_Flt[10] = interfaceHeight;

   AuxArray_Int[ 0] = Jet_Ambient;

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
//                GenePtr           : Array of pointers for general potential tables
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_IsothermalSlab( const double x, const double y, const double z, const double Time,
                                   const double UserArray_Flt[], const int UserArray_Int[],
                                   const ExtPotUsage_t Usage, const real PotTable[], void **GenePtr )
{

// halo potential
   const double cz                                =       UserArray_Flt[ 2];   // z ...
   const real   IsothermalSlab_VelocityDispersion = (real)UserArray_Flt[ 3];   //
   const real   IsothermalSlab_PeakDens           = (real)UserArray_Flt[ 4];   //
   const real   NewtonG                           =       UserArray_Flt[ 5];
   const real   BoxSize_Z                         =       UserArray_Flt[ 6];
   const real   IsothermalSlab_Truncation         = (real)UserArray_Flt[ 7];   //
   const real   distance_h                        = (real)UserArray_Flt[ 8];
   const real   v_halo                            = (real)UserArray_Flt[ 9];
   const real   interfaceHeight                   = (real)UserArray_Flt[10];
   const int    Jet_Ambient                       =       UserArray_Int[ 0];

   const double IsothermalSlab_VelocityDispersion_Sqr = SQR(IsothermalSlab_VelocityDispersion);
   const double dz = z - cz;

   real stellarDiskPot = sqrt( ( 2.0*M_PI*NewtonG*IsothermalSlab_PeakDens ) / IsothermalSlab_VelocityDispersion_Sqr );

#  ifndef __CUDACC__
   if ( cz != 0.5*BoxSize_Z )   Aux_Error( ERROR_INFO, "We expect the z-position of stellar disk is at box-center!\n" );
#  endif

   real LogPot;

   if ( Jet_Ambient == 2 )
   {
      if ( fabs(dz) > IsothermalSlab_Truncation )
      {
         stellarDiskPot = 2.0*IsothermalSlab_VelocityDispersion_Sqr*log(cosh(IsothermalSlab_Truncation*stellarDiskPot));
         LogPot = SQR(v_halo) * log(SQR(IsothermalSlab_Truncation) + SQR(distance_h));
      } else
      {
         stellarDiskPot = 2.0*IsothermalSlab_VelocityDispersion_Sqr*log(cosh(dz*stellarDiskPot));
         LogPot = SQR(v_halo) * log(dz*dz + SQR(distance_h));
      }
   }
   else if ( Jet_Ambient == 3 )
   {
      if ( fabs(dz) > interfaceHeight )
      {
         stellarDiskPot = 2.0*IsothermalSlab_VelocityDispersion_Sqr*log(cosh(interfaceHeight*stellarDiskPot));
         LogPot = SQR(v_halo) * log(SQR(interfaceHeight) + SQR(distance_h));
      } else
      {
         stellarDiskPot = 2.0*IsothermalSlab_VelocityDispersion_Sqr*log(cosh(dz*stellarDiskPot));
         LogPot = SQR(v_halo) * log(dz*dz + SQR(distance_h));
      }
   } // if ( Jet_Ambient == 2 ) ... else if ... else ...

   real TotPot = stellarDiskPot + LogPot;

   return TotPot;

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
