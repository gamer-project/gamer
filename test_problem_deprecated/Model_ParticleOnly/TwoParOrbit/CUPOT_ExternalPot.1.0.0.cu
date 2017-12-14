#ifdef __CUDACC__
#include "Macro.h"
#else
#include "GAMER.h"
#endif
#include "CUPOT.h"

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  CUPOT_ExternalPot / CPU_ExternalPot
// Description :  1. Cacalculate the external potential from the input coordinates and time
//                2. This function will be invoked by both the CPU and GPU codes
//
// Parameter   :  x/y/z     : Spatial coordinates
//                Time      : Current physical time
//                UserArray : User-provided auxiliary array (set by "Init_ExternalPot")
//
// Return      :  external potential
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__device__
real CUPOT_ExternalPot( const double x, const double y, const double z, const double Time, const double UserArray[] )
#else
real   CPU_ExternalPot( const double x, const double y, const double z, const double Time, const double UserArray[] )
#endif
{

   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
   const real   GM_4   = (real)UserArray[3];   // 0.25*G*m
   const real   dx     = (real)(x - Cen[0]);
   const real   dy     = (real)(y - Cen[1]);
   const real   dz     = (real)(z - Cen[2]);
   const real   _r     = (real)1.0/SQRT( dx*dx + dy*dy + dz*dz );

   return -GM_4*_r;

} // FUNCTION : CUPOT_ExternalPot // CPU_ExternalPot



#endif // #ifdef GRAVITY
