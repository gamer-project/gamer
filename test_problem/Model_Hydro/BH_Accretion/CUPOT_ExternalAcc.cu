#include "Copyright.h"
#ifdef __CUDACC__
#include "Macro.h"
#else
#include "GAMER.h"
#endif
#include "CUPOT.h"

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  CUPOT_ExternalAcc / CPU_ExternlAcc
// Description :  1. Cacalculate the external acceleration from the input coordinates and time
//                2. This function will be invoked in both the CPU and GPU codes
//                3. "__forceinline__" is required since this device function will be invoked by more than one kernel
//                   (e.g., CUPOT_HydroGravitySolver, CUFLU_ComputeFlux <-- which will be called by different fluid solvers)
//
// Parameter   :  Acc         : Array to store the output external acceleration
//                x/y/z       : Spatial coordinates
//                Time        : Current physical time
//                UserArray   : User-provided auxiliary array (set by "Init_ExternalPot")
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

   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
   const real GM       = (real)UserArray[3];
   const real eps      = (real)UserArray[4];
   const real dx       = (real)(x - Cen[0]);
   const real dy       = (real)(y - Cen[1]);
   const real dz       = (real)(z - Cen[2]);
   const real r        = SQRT( dx*dx + dy*dy + dz*dz );

//#  define SOFTEN_PLUMMER
#  define SOFTEN_RUFFERT

// Plummer
#  if   defined SOFTEN_PLUMMER
   const real _r3      = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps), (real)-1.5 );

// Ruffert 1994
#  elif defined SOFTEN_RUFFERT
   const real tmp      = EXP( -SQR(r)/SQR(eps) );
   const real _r3      = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps)*tmp, (real)-1.5 )*( (real)1.0 - tmp );

#  else
#  error : UNKNOWN SOFTEN !!
#  endif

   Acc[0] = -GM*_r3*dx;
   Acc[1] = -GM*_r3*dy;
   Acc[2] = -GM*_r3*dz;

} // FUNCTION : CUPOT_ExternalAcc / CPU_ExternalAcc



#endif // #ifdef GRAVITY
