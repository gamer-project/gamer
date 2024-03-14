#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#endif

#ifdef GRAVITY




// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__
extern double   HaloMerger_ExtPot_UniDenSph_CenCoordX;
extern double   HaloMerger_ExtPot_UniDenSph_CenCoordY;
extern double   HaloMerger_ExtPot_UniDenSph_CenCoordZ;
extern double   HaloMerger_ExtPot_UniDenSph_VelocityX;
extern double   HaloMerger_ExtPot_UniDenSph_VelocityY;
extern double   HaloMerger_ExtPot_UniDenSph_VelocityZ;
extern double   HaloMerger_ExtPot_UniDenSph_M;
extern double   HaloMerger_ExtPot_UniDenSph_R;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_ELBDM_HaloMerger
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_ELBDM_HaloMerger()
//
// Note        :  1. Invoked by Init_ExtPot_ELBDM_HaloMerger()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//                Time             : Target physical time
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_ELBDM_HaloMerger( double AuxArray_Flt[], int AuxArray_Int[], const double Time )
{

   AuxArray_Flt[0] = HaloMerger_ExtPot_UniDenSph_CenCoordX;
   AuxArray_Flt[1] = HaloMerger_ExtPot_UniDenSph_CenCoordY;
   AuxArray_Flt[2] = HaloMerger_ExtPot_UniDenSph_CenCoordZ;
   AuxArray_Flt[3] = HaloMerger_ExtPot_UniDenSph_VelocityX;
   AuxArray_Flt[4] = HaloMerger_ExtPot_UniDenSph_VelocityY;
   AuxArray_Flt[5] = HaloMerger_ExtPot_UniDenSph_VelocityZ;
   AuxArray_Flt[6] = NEWTON_G*HaloMerger_ExtPot_UniDenSph_M;
   AuxArray_Flt[7] = HaloMerger_ExtPot_UniDenSph_R;

} // FUNCTION : SetExtPotAuxArray_ELBDM_HaloMerger
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_ELBDM_HaloMerger
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_ELBDM_HaloMerger(), where
//                      UserArray_Flt[0] = x coordinate of the external potential center
//                      UserArray_Flt[1] = y ...
//                      UserArray_Flt[2] = z ..
//                      UserArray_Flt[3] = x component of the velocity of the external potential center
//                      UserArray_Flt[4] = y component of the velocity of the external potential center
//                      UserArray_Flt[5] = z component of the velocity of the external potential center
//                      UserArray_Flt[6] = gravitational_constant*sphere_mass
//                      UserArray_Flt[7] = sphere_radius
//                3. Currently it does not support the soften length
//                4. GenePtr has the size of EXT_POT_NGENE_MAX defined in Macro.h (default = 6)
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
static real ExtPot_ELBDM_HaloMerger( const double x, const double y, const double z, const double Time,
                                     const double UserArray_Flt[], const int UserArray_Int[],
                                     const ExtPotUsage_t Usage, const real PotTable[], void **GenePtr )
{

   const double Cen[3] = { UserArray_Flt[0], UserArray_Flt[1], UserArray_Flt[2] };
   const double Vel[3] = { UserArray_Flt[3], UserArray_Flt[4], UserArray_Flt[5] };
   const real   GM     = (real)UserArray_Flt[6];
   const real   R      = (real)UserArray_Flt[7];
   const real   dx     = (real)(x - Cen[0] - Vel[0]*Time);
   const real   dy     = (real)(y - Cen[1] - Vel[1]*Time);
   const real   dz     = (real)(z - Cen[2] - Vel[2]*Time);
   const real   r      = SQRT( dx*dx + dy*dy + dz*dz );
   const real   _r     = (real)1.0/r;
   const real   GM_2R3 = GM/((real)2.0*CUBE(R));
   const real   Offset = (real)0.75*(GM/R);  // Shift the potential to reduce the absolute value of potential
                                             // and thus increase the time step

   if ( r >= R )
      return -GM*_r + Offset;
   else
      return GM_2R3*(r*r-(real)3.0*R*R) + Offset;

} // FUNCTION : ExtPot_ELBDM_HaloMerger



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_ELBDM_HaloMerger;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_ELBDM_HaloMerger
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_ELBDM_HaloMerger()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtPot_ELBDM_HaloMerger( ExtPot_t &CPUExtPot_Ptr, ExtPot_t &GPUExtPot_Ptr )
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_ELBDM_HaloMerger( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_ELBDM_HaloMerger( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_ELBDM_HaloMerger( double [], int [], const double );
void SetCPUExtPot_ELBDM_HaloMerger( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_ELBDM_HaloMerger( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_ELBDM_HaloMerger
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
void Init_ExtPot_ELBDM_HaloMerger()
{

   SetExtPotAuxArray_ELBDM_HaloMerger( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int, Time[0] );
   SetCPUExtPot_ELBDM_HaloMerger( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_ELBDM_HaloMerger( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_ELBDM_HaloMerger

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
