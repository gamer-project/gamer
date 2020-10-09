#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY




// =================================
// I. Set an auxiliary array
// =================================

#ifndef __CUDACC__

extern double Plummer_Rho0;
extern double Plummer_R0;
extern double Plummer_Center[3];
extern double Plummer_ExtPotMFrac;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_Plummer
// Description :  Set the auxiliary array ExtPot_AuxArray[] used by ExtPot_Plummer()
//
// Note        :  1. Invoked by Init_ExtPot_Plummer()
//                2. AuxArray[] has the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_Plummer( double AuxArray[] )
{

// potential = -G*M/(r^2+R0^2)^(1/2)
   const double ExtPotM = (4.0/3.0)*M_PI*CUBE(Plummer_R0)*Plummer_Rho0*Plummer_ExtPotMFrac;

   AuxArray[0] = Plummer_Center[0];    // x coordinate of the Plummer center
   AuxArray[1] = Plummer_Center[1];    // y ...
   AuxArray[2] = Plummer_Center[2];    // z ...
   AuxArray[3] = SQR( Plummer_R0 );    // scale_radius^2
   AuxArray[4] = -NEWTON_G*ExtPotM;    // -G*M

} // FUNCTION : SetExtPotAuxArray_Plummer
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_Plummer
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by SetExtPotAuxArray_Plummer()
//
// Parameter   :  x/y/z     : Target spatial coordinates
//                Time      : Target physical time
//                UserArray : User-provided auxiliary array
//                Usage     : Different usages of external potential when computing total potential on level Lv
//                            --> EXT_POT_USAGE_ADD     : add external potential on Lv
//                                EXT_POT_USAGE_SUB     : subtract external potential for preparing self-gravity potential on Lv-1
//                                EXT_POT_USAGE_SUB_TINT: like SUB but for temporal interpolation
//                            --> This parameter is useless in most cases
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_Plummer( const double x, const double y, const double z, const double Time, const double UserArray[],
                            const ExtPotUsage_t Usage )
{

// potential = -G*M/(r^2+R0^2)^(1/2)
   const double cx  =       UserArray[0];    // x coordinate of the Plummer center
   const double cy  =       UserArray[1];    // y ...
   const double cz  =       UserArray[2];    // z ...
   const real   a2  = (real)UserArray[3];    // scale_radius^2
   const real   mGM = (real)UserArray[4];    // -G*M

   const real   dx  = (real)(x - cx);
   const real   dy  = (real)(y - cy);
   const real   dz  = (real)(z - cz);
   const real   r2  = SQR(dx) + SQR(dy) + SQR(dz);
   const real   pot = mGM / SQRT( r2 + a2 );

   return pot;

} // FUNCTION : ExtPot_Plummer



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_Plummer;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_Plummer
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_Plummer()
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_Plummer( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_Plummer( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_Plummer( double [] );
void SetCPUExtPot_Plummer( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_Plummer( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_Plummer
// Description :  Initialize external potential
//
// Note        :  1. Set an auxiliary array by invoking SetExtPotAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external potential major routines by invoking SetCPU/GPUExtPot_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtPot_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtPot_Plummer()
{

   SetExtPotAuxArray_Plummer( ExtPot_AuxArray );
   SetCPUExtPot_Plummer( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_Plummer( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_Plummer

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
