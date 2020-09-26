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
extern double Plummer_ExtMFrac;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtAccAuxArray_Plummer
// Description :  Set the auxiliary array ExtAcc_AuxArray[] used by ExtAcc_Plummer()
//
// Note        :  1. Invoked by Init_ExtAcc_Plummer()
//                2. AuxArray[] has the size of EXT_ACC_NAUX_MAX defined in Macro.h (default = 10)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void SetExtAccAuxArray_Plummer( double AuxArray[] )
{

// acceleration = -G*Mtot*r/(r^2+R0^2)^(3/2)
   const double Mtot = (4.0/3.0)*M_PI*CUBE(Plummer_R0)*Plummer_Rho0*Plummer_ExtMFrac;

   AuxArray[0] = Plummer_Center[0];    // x coordinate of the external acceleration center
   AuxArray[1] = Plummer_Center[1];    // y ...
   AuxArray[2] = Plummer_Center[2];    // z ...
   AuxArray[3] = SQR( Plummer_R0 );    // scale_radius^2
   AuxArray[4] = -NEWTON_G*Mtot;       // -G*total_mass

} // FUNCTION : SetExtAccAuxArray_Plummer
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external acceleration
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_Plummer
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by SetExtAccAuxArray_Plummer(), where
//                      UserArray[0] = x coordinate of the external acceleration center
//                      UserArray[1] = y ...
//                      UserArray[2] = z ..
//                      UserArray[3] = gravitational_constant*point_source_mass
//                      UserArray[4] = soften_length (<=0.0 --> disable)
//                3. Two different soften length implementations are supported
//                   --> SOFTEN_PLUMMER & SOFTEN_RUFFERT
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Target physical time
//                UserArray : User-provided auxiliary array
//
// Return      :  External acceleration Acc[] at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void ExtAcc_Plummer( real Acc[], const double x, const double y, const double z, const double Time,
                              const double UserArray[] )
{

// acceleration = -G*Mtot*r/(r^2+a^2)^(3/2)
   const double cx     =       UserArray[0];    // x coordinate of the external acceleration center
   const double cy     =       UserArray[1];    // y ...
   const double cz     =       UserArray[2];    // z ...
   const real   a2     = (real)UserArray[3];    // scale_radius^2
   const real   mGMtot = (real)UserArray[4];    // -G*total_mass

   const real   dx     = (real)(x - cx);
   const real   dy     = (real)(y - cy);
   const real   dz     = (real)(z - cz);
   const real   r2     = SQR(dx) + SQR(dy) + SQR(dz);
   const real   tmp    = mGMtot * POW( r2+a2, (real)-1.5 );

   Acc[0] = tmp*dx;
   Acc[1] = tmp*dy;
   Acc[2] = tmp*dz;

} // FUNCTION : ExtAcc_Plummer



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtAcc_t ExtAcc_Ptr = ExtAcc_Plummer;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_Plummer
// Description :  Return the function pointers of the CPU/GPU external acceleration routines
//
// Note        :  1. Invoked by Init_ExtAcc_Plummer()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtAcc_Plummer( ExtAcc_t &CPUExtAcc_Ptr, ExtAcc_t &GPUExtAcc_Ptr )
//
// Parameter   :  CPU/GPUExtAcc_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtAcc_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtAcc_Plummer( ExtAcc_t &GPUExtAcc_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtAcc_Ptr, ExtAcc_Ptr, sizeof(ExtAcc_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtAcc_Plummer( ExtAcc_t &CPUExtAcc_Ptr )
{
   CPUExtAcc_Ptr = ExtAcc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtAccAuxArray_Plummer( double [] );
void SetCPUExtAcc_Plummer( ExtAcc_t & );
#ifdef GPU
void SetGPUExtAcc_Plummer( ExtAcc_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtAcc_Plummer
// Description :  Initialize external acceleration
//
// Note        :  1. Set an auxiliary array by invoking SetExtAccAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external acceleration major routines by invoking SetCPU/GPUExtAcc_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtAcc_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  OnlySetAuxArray : See "src/SelfGravity/Init_ExtAccPot.cpp"
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtAcc_Plummer( const bool OnlySetAuxArray )
{

   SetExtAccAuxArray_Plummer( ExtAcc_AuxArray );

   if ( ! OnlySetAuxArray )
   {
      SetCPUExtAcc_Plummer( CPUExtAcc_Ptr );
#     ifdef GPU
      SetGPUExtAcc_Plummer( GPUExtAcc_Ptr );
#     endif
   }

} // FUNCTION : Init_ExtAcc_Plummer

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
