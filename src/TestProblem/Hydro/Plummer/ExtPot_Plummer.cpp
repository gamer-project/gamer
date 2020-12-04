#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY




// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__

extern double Plummer_Rho0;
extern double Plummer_R0;
extern double Plummer_Center[3];
extern double Plummer_ExtPotMFrac;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_Plummer
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_Plummer()
//
// Note        :  1. Invoked by Init_ExtPot_Plummer()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_Plummer( double AuxArray_Flt[], int AuxArray_Int[] )
{

// potential = -G*M/(r^2+R0^2)^(1/2)
   const double ExtPotM = (4.0/3.0)*M_PI*CUBE(Plummer_R0)*Plummer_Rho0*Plummer_ExtPotMFrac;

   AuxArray_Flt[0] = Plummer_Center[0];   // x coordinate of the Plummer center
   AuxArray_Flt[1] = Plummer_Center[1];   // y ...
   AuxArray_Flt[2] = Plummer_Center[2];   // z ...
   AuxArray_Flt[3] = SQR( Plummer_R0 );   // scale_radius^2
   AuxArray_Flt[4] = -NEWTON_G*ExtPotM;   // -G*M

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
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_Plummer()
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
static real ExtPot_Plummer( const double x, const double y, const double z, const double Time,
                            const double UserArray_Flt[], const int UserArray_Int[],
                            const ExtPotUsage_t Usage, const real PotTable[] )
{

// potential = -G*M/(r^2+R0^2)^(1/2)
   const double cx  =       UserArray_Flt[0];   // x coordinate of the Plummer center
   const double cy  =       UserArray_Flt[1];   // y ...
   const double cz  =       UserArray_Flt[2];   // z ...
   const real   a2  = (real)UserArray_Flt[3];   // scale_radius^2
   const real   mGM = (real)UserArray_Flt[4];   // -G*M

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
void SetExtPotAuxArray_Plummer( double [], int [] );
void SetCPUExtPot_Plummer( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_Plummer( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_Plummer
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
void Init_ExtPot_Plummer()
{

   SetExtPotAuxArray_Plummer( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int );
   SetCPUExtPot_Plummer( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_Plummer( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_Plummer

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
