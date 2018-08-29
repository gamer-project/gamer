#ifdef __CUDACC__
#include "Macro.h"
#else
#include "GAMER.h"
#endif
#include "CUPOT.h"

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  CUPOT_ExternalAcc / CPU_ExternlAcc
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function will be invoked by both CPU and GPU
//                2. "__forceinline__" is required since this device function will be invoked
//                   by more than one kernels (e.g., CUPOT_HydroGravitySolver, CUFLU_ComputeFlux)
//                3. The auxiliary array "UserArray" is set by "Init_ExternalAcc_Ptr", which
//                   points to "Init_ExternalAcc()" by default but may be overwritten by various
//                   test problem initializers
//                4. By default we assume
//                     UserArray[0] = x coordinate of the external acceleration center
//                     UserArray[1] = y ...
//                     UserArray[2] = z ..
//                     UserArray[3] = gravitational_constant*point_source_mass
//                     UserArray[4] = soften_length (<=0.0 --> disable)
//                   --> but one can easily modify this file to change the default behavior
//                5. Two different soften length implementations are supported
//                   --> SOFTEN_PLUMMER & SOFTEN_RUFFERT
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Current physical time
//                UserArray : User-provided auxiliary array (set by "Init_ExternalAcc_Ptr")
//
// Return      :  Acc
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__forceinline__ __device__
void CUPOT_ExternalAcc( real Acc[], const double x, const double y, const double z, const double Time, const double UserArray[] )
#else
void   CPU_ExternalAcc( real Acc[], const double x, const double y, const double z, const double Time, const double UserArray[] )
#endif
{

   const double Cen[3]  = { UserArray[0], UserArray[1], UserArray[2] };    // cluster center
   const real   coeff   = (real)UserArray[3];                              // -G*M200/( log(1+c) - c/(1+c) )
   const real   r0      = (real)UserArray[4];                              // scale radius
   const real   Vy      = (real)UserArray[5];                              // jet infall velocity along y
   const real   dx      = (real)0.0;                                       // do not consider dx for now
   const real   dy      = (real)(y - Cen[1] - Vy*Time);
   const real   dz      = (real)0.0;                                       // do not consider dz for now
   const real   r       = SQRT( dx*dx + dy*dy + dz*dz );
   const real   s       = r/r0;
   const real   one     = (real)1.0;
   const real   _r3     = one/CUBE(r);
   const real   acc_r   = coeff*( LOG(one+s) - s/(one+s) )*_r3;

   Acc[0] = acc_r*dx;
   Acc[1] = acc_r*dy;
   Acc[2] = acc_r*dz;

} // FUNCTION : CUPOT_ExternalAcc / CPU_ExternalAcc



#endif // #ifdef GRAVITY
