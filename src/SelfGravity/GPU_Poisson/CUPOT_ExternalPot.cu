#ifdef __CUDACC__
#include "Macro.h"
#else
#include "GAMER.h"
#endif
#include "CUPOT.h"

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  CUPOT_ExternalPot / CPU_ExternalPot
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function will be invoked by both CPU and GPU
//                2. The auxiliary array "UserArray" is set by "Init_ExternalPot_Ptr", which
//                   points to "Init_ExternalPot()" by default but may be overwritten by various
//                   test problem initializers
//                3. By default we assume
//                     UserArray[0] = x coordinate of the external acceleration center
//                     UserArray[1] = y ...
//                     UserArray[2] = z ..
//                     UserArray[3] = gravitational_constant*point_source_mass
//                   --> but one can easily modify this file to change the default behavior
//                4. Currently it does not support the soften length
//
// Return      :  External potential
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__device__
real CUPOT_ExternalPot( const double x, const double y, const double z, const double Time, const double UserArray[] )
#else
real   CPU_ExternalPot( const double x, const double y, const double z, const double Time, const double UserArray[] )
#endif
{

   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
   const real   GM     = (real)UserArray[3];
   const real   dx     = (real)(x - Cen[0]);
   const real   dy     = (real)(y - Cen[1]);
   const real   dz     = (real)(z - Cen[2]);
   const real   _r     = 1.0/SQRT( dx*dx + dy*dy + dz*dz );

   return -GM*_r;

} // FUNCTION : CUPOT_ExternalPot // CPU_ExternalPot



#endif // #ifdef GRAVITY
