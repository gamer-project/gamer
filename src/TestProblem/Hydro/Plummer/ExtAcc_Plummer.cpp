#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_Plummer
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Current physical time
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
// get the CPU/GPU function pointers
// =================================

#ifdef __CUDACC__
__device__
#endif
static ExtAcc_t ExtAcc_Ptr = ExtAcc_Plummer;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_Plummer
// Description :  Return the function pointers to the CPU/GPU external acceleration routines
//
// Note        :  None
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



#endif // #ifdef GRAVITY
