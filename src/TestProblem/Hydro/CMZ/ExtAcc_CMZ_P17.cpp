#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_BarredPot
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
static void ExtAcc_BarredPot( real Acc[], const double x, const double y, const double z, const double Time,
                            const double UserArray[] )
{

   const double cx     =       UserArray[5];    // x coordinate of the external acceleration center
   const double cy     =       UserArray[6];    // y ...
   const double cz     =       UserArray[7];    // z ...

   const double V02    =       UserArray[0];
   const double q2     =       UserArray[1];
   const double Rc2    =       UserArray[2];
   const double Omegab =       UserArray[3];
   const double fullBS =       UserArray[4];

// convert x,y,z to be centered on mesh center

   const real   dx     = (real)(x - cx);
   const real   dy     = (real)(y - cy);
   const real   dz     = (real)(z - cz);

   const double angle  =  -1.*Time*Omegab;
   const double dxrot  =  dx*cos(angle) - dy*sin(angle);
   const double dyrot  =  dx*sin(angle) + dy*cos(angle);

   const real   A      = (real)(Rc2 + SQR(dxrot) +(SQR(dyrot) + SQR(dz))/q2);
   const real   B      = (real)(Rc2 + SQR(dxrot) +(SQR(dyrot) + SQR(dz)));

// const real   frac   = std::min(Time/fullBS,1.0); 
   const real   frac   = FMIN(Time/fullBS,1.0); 
 
   const real fxs = real(-V02*(dx/B));
   const real fys = real(-V02*(dy/B));
   const real fzs = real(-V02*(dz/B));
     
   const real fxr = real(-V02*(dxrot/A));
   const real fyr = real(-V02*(dyrot/A)/q2);

   const real fxb = real(fxr*cos(-angle) - fyr*sin(-angle));
   const real fyb = real(fxr*sin(-angle) + fyr*cos(-angle));
   const real fzb = real(-V02*(dz/A   )/q2);


   Acc[0] = fxb*frac + fxs*(1-frac);
   Acc[1] = fyb*frac + fys*(1-frac);
   Acc[2] = fzb*frac + fzs*(1-frac);

} // FUNCTION : ExtAcc_BarredPot



// =================================
// get the CPU/GPU function pointers
// =================================

#ifdef __CUDACC__
__device__
#endif
static ExtAcc_t ExtAcc_Ptr = ExtAcc_BarredPot;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_BarredPot
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
void SetGPUExtAcc_BarredPot( ExtAcc_t &GPUExtAcc_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtAcc_Ptr, ExtAcc_Ptr, sizeof(ExtAcc_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtAcc_BarredPot( ExtAcc_t &CPUExtAcc_Ptr )
{
   CPUExtAcc_Ptr = ExtAcc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

extern double BarredPot_V0;
extern double BarredPot_q;
extern double BarredPot_Rc;
extern double BarredPot_Omegabar;
extern double BarredPot_fullBS;
extern double BarredPot_MeshCenter[3];

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtAccAuxArray_BarredPot
// Description :  Set the auxiliary array ExtAcc_AuxArray[] used by ExtAcc_BarredPot()
//
// Note        :  1. AuxArray[] has the size of EXT_ACC_NAUX_MAX defined in Macro.h (default = 10)
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void Init_ExtAccAuxArray_BarredPot( double AuxArray[] )
{

   AuxArray[0] = SQR(BarredPot_V0);    // amplitude^2  
   AuxArray[1] = SQR(BarredPot_q );    // axis ratio^2
   AuxArray[2] = SQR(BarredPot_Rc);    // core radius^2
   AuxArray[3] = BarredPot_Omegabar;   // bar pattern speed
   AuxArray[4] = BarredPot_fullBS;     // bar growth time

   AuxArray[5] = BarredPot_MeshCenter[0];    // x coordinate of the external acceleration center
   AuxArray[6] = BarredPot_MeshCenter[1];    // y ...
   AuxArray[7] = BarredPot_MeshCenter[2];    // z ...



} // FUNCTION : Init_ExtAccAuxArray_BarredPot
#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
