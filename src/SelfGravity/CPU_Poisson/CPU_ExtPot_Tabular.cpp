#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY



/********************************************************
1. Tabular external potential

2. This file is shared by both CPU and GPU

   GPU_Poisson/CUPOT_ExtPot_Tabular.cu -> CPU_Poisson/CPU_ExtPot_Tabular.cpp

3. Three steps are required to implement external potential

   I.   Set auxiliary arrays
        --> SetExtPotAuxArray_Tabular()

   II.  Specify external potential
        --> ExtPot_Tabular()

   III. Set initialization functions
        --> SetGPUExtPot_Tabular()
            SetCPUExtPot_Tabular()
            Init_ExtPot_Tabular()

4. The external potential major routine, ExtPot_Tabular(),
   must be thread-safe and not use any global variable

5. Reference: https://github.com/gamer-project/gamer/wiki/Gravity#external-accelerationpotential
********************************************************/



// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_Tabular
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_Tabular()
//
// Note        :  1. Invoked by Init_ExtPot_Tabular()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_Tabular( double AuxArray_Flt[], int AuxArray_Int[] )
{

// floating-point parameters
   AuxArray_Flt[0] = EXT_POT_TABLE_EDGEL[0];    // left-edge x/y/z coordinates of the table
   AuxArray_Flt[1] = EXT_POT_TABLE_EDGEL[1];
   AuxArray_Flt[2] = EXT_POT_TABLE_EDGEL[2];
   AuxArray_Flt[3] = 1.0 / EXT_POT_TABLE_DH;    // 1/dh

// integer parameters
   AuxArray_Int[0] = EXT_POT_TABLE_NPOINT[0];   // table sizes along x/y/z
   AuxArray_Int[1] = EXT_POT_TABLE_NPOINT[1];
   AuxArray_Int[2] = EXT_POT_TABLE_NPOINT[2];

} // FUNCTION : SetExtPotAuxArray_Tabular
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_Tabular
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_Tabular()
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
static real ExtPot_Tabular( const double x, const double y, const double z, const double Time,
                            const double UserArray_Flt[], const int UserArray_Int[],
                            const ExtPotUsage_t Usage, const real PotTable[] )
{

   const double EdgeL_x  = UserArray_Flt[0];
   const double EdgeL_y  = UserArray_Flt[1];
   const double EdgeL_z  = UserArray_Flt[2];
   const double _dh      = UserArray_Flt[3];

   const int    NPoint_x = UserArray_Int[0];
   const int    NPoint_y = UserArray_Int[1];
   const int    didx_x   = 1;
   const int    didx_y   = NPoint_x;
   const int    didx_z   = NPoint_x*NPoint_y;

   const real   ONE      = (real)1.0;


// 1. get the lower corner cell index
   int  idx_x, idx_y, idx_z;
   long idx0;
   real dx, dy, dz;

   dx    = real( ( x - EdgeL_x )*_dh );
   dy    = real( ( y - EdgeL_y )*_dh );
   dz    = real( ( z - EdgeL_z )*_dh );
   idx_x = int( dx );
   idx_y = int( dy );
   idx_z = int( dz );
   idx0  = long( idx_x*didx_x + idx_y*didx_y ) + (long)idx_z*didx_z;

// it should never happen even considering round-off errors!
#  if ( defined GAMER_DEBUG  &&  !defined __CUDACC__ )
   const int NPoint_z = UserArray_Int[2];

   if ( idx_x < 0  ||  idx_x+1 >= NPoint_x )
      Aux_Error( ERROR_INFO, "x index outside the table range (x %14.7e, EdgeL %14.7e, _dh %13.7e, idx %d) !!\n",
                 x, EdgeL_x, _dh, idx_x );

   if ( idx_y < 0  ||  idx_y+1 >= NPoint_y )
      Aux_Error( ERROR_INFO, "y index outside the table range (y %14.7e, EdgeL %14.7e, _dh %13.7e, idx %d) !!\n",
                 y, EdgeL_y, _dh, idx_y );

   if ( idx_z < 0  ||  idx_z+1 >= NPoint_z )
      Aux_Error( ERROR_INFO, "z index outside the table range (z %14.7e, EdgeL %14.7e, _dh %13.7e, idx %d) !!\n",
                 z, EdgeL_z, _dh, idx_z );
#  endif


// 2. compute potential by trilinear interpolation
   real weight_xL, weight_yL, weight_zL;
   real weight_xR, weight_yR, weight_zR;
   real ExtPot;

   weight_xR = dx - (real)idx_x;
   weight_yR = dy - (real)idx_y;
   weight_zR = dz - (real)idx_z;
   weight_xL = ONE - weight_xR;
   weight_yL = ONE - weight_yR;
   weight_zL = ONE - weight_zR;

   ExtPot = PotTable[ idx0                            ] * weight_xL * weight_yL * weight_zL +
            PotTable[ idx0 + didx_x                   ] * weight_xR * weight_yL * weight_zL +
            PotTable[ idx0          + didx_y          ] * weight_xL * weight_yR * weight_zL +
            PotTable[ idx0                   + didx_z ] * weight_xL * weight_yL * weight_zR +
            PotTable[ idx0 + didx_x + didx_y          ] * weight_xR * weight_yR * weight_zL +
            PotTable[ idx0          + didx_y + didx_z ] * weight_xL * weight_yR * weight_zR +
            PotTable[ idx0 + didx_x          + didx_z ] * weight_xR * weight_yL * weight_zR +
            PotTable[ idx0 + didx_x + didx_y + didx_z ] * weight_xR * weight_yR * weight_zR;

   return ExtPot;

} // FUNCTION : ExtPot_Tabular



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_Tabular;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_Tabular
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_Tabular()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtPot_Tabular( ExtPot_t &CPUExtPot_Ptr, ExtPot_t &GPUExtPot_Ptr )
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_Tabular( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_Tabular( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_Tabular( double [], int [] );
void SetCPUExtPot_Tabular( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_Tabular( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_Tabular
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
void Init_ExtPot_Tabular()
{

   SetExtPotAuxArray_Tabular( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int );
   SetCPUExtPot_Tabular( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_Tabular( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_Tabular

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
