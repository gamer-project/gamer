#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#endif

#ifdef GRAVITY

// =================================
// // I. Set auxiliary arrays
// // =================================
//
#ifndef __CUDACC__

extern double BarredPot_V0;
extern double BarredPot_q;
extern double BarredPot_Rc;
extern double BarredPot_Omegabar;
extern double BarredPot_fullBS;
extern double BarredPot_MeshCenter[3];

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_BarredPot
// Description :  Set the auxiliary array ExtPot_AuxArray_Flt/Int[] used by ExtPot_BarredPot()
//
// Note        :  1. Invoked by Init_ExtPot_BarredPot()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//                Time             : Target physical time
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_BarredPot(double AuxArray_Flt[], int AuxArray_Int[], const double Time )
{

   AuxArray_Flt[0] = SQR(BarredPot_V0);    // amplitude^2
   AuxArray_Flt[1] = SQR(BarredPot_q );    // axis ratio^2
   AuxArray_Flt[2] = SQR(BarredPot_Rc);    // core radius^2
   AuxArray_Flt[3] = BarredPot_Omegabar;   // bar pattern speed
   AuxArray_Flt[4] = BarredPot_fullBS;     // bar growth time

   AuxArray_Flt[5] = BarredPot_MeshCenter[0];    // x coordinate of the external acceleration center
   AuxArray_Flt[6] = BarredPot_MeshCenter[1];    // y ...
   AuxArray_Flt[7] = BarredPot_MeshCenter[2];    // z ...



} // FUNCTION : SetExtPotAuxArray_BarredPot
#endif // #ifndef __CUDACC__




// =================================
// II. Specify external potential
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_BarredPot
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_BarredPot()
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
static real ExtPot_BarredPot(const double x, const double y, const double z, const double Time,
                             const double UserArray_Flt[], const int UserArray_Int[],
                             const ExtPotUsage_t Usage, const real PotTable[], void **GenePtr )
{

   const double cx     =       UserArray_Flt[5];    // x coordinate of the external acceleration center
   const double cy     =       UserArray_Flt[6];    // y ...
   const double cz     =       UserArray_Flt[7];    // z ...

   const double V02    =       UserArray_Flt[0];
   const double q2     =       UserArray_Flt[1];
   const double Rc2    =       UserArray_Flt[2];
   const double Omegab =       UserArray_Flt[3];
   const double fullBS =       UserArray_Flt[4];

// convert x,y,z to be centered on mesh center

   const real   dx     = (real)(x - cx);
   const real   dy     = (real)(y - cy);
   const real   dz     = (real)(z - cz);
   const real   q      = (real)sqrt(q2);

   const double angle  =  -1.*Time*Omegab;
   const double dxrot  =  dx*cos(angle) - dy*sin(angle);
   const double dyrot  =  dx*sin(angle) + dy*cos(angle);

   const real   A      = (real)(Rc2 + SQR(dxrot) +(SQR(dyrot) + SQR(dz))/q2);
   const real   B      = (real)(Rc2 + SQR(dxrot) +(SQR(dyrot) + SQR(dz)));
   const real   pconst = (real)2.*log((1+q)*(1-q)/(4*q));

// const real   frac   = std::min(Time/fullBS,1.0);
   const real   frac   = FMIN(Time/fullBS,1.0);
   const real   pots   = (real)0.5*V02*(log(B)-pconst);
   const real   potb   = (real)0.5*V02*(log(A)-pconst);

   const real   pot    = potb*frac + pots*(1-frac);

   return pot;

} // FUNCTION : ExtPot_BarredPot



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_BarredPot;


//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_BarredPot
// Description :  Return the function pointers of the CPU/GPU external acceleration routines
//
// Note        :  1. Invoked by Init_ExtPot_BarredPot()
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------


#ifdef __CUDACC__
__host__

void SetGPUExtPot_BarredPot( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_BarredPot( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__


// local function prototypes
void SetExtPotAuxArray_BarredPot( double [], int [], const double );
void SetCPUExtPot_BarredPot( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_BarredPot( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_BarredPot
// Description :  Initialize external acceleration
//
// Note        :  1. Set an auxiliary array by invoking SetExtPotAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external acceleration major routines by invoking SetCPU/GPUExtPot_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtPot_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtPot_BarredPot()
{

   SetExtPotAuxArray_BarredPot( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int, Time[0] );

   SetCPUExtPot_BarredPot( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_BarredPot( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_BarredPot

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
