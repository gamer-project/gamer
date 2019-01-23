#ifndef __EXTERNALPOT__
#define __EXTERNALPOT__



#include "CUPOT.h"

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  ExternalPot
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by "Init_ExternalPot_Ptr", which
//                   points to Init_ExternalPot() by default but may be overwritten by various
//                   test problem initializers
//                3. By default we assume
//                     UserArray[0] = x coordinate of the external acceleration center
//                     UserArray[1] = y ...
//                     UserArray[2] = z ..
//                     UserArray[3] = gravitational_constant*point_source_mass
//                   --> But one can easily modify this file to change the default behavior
//                4. Currently it does not support the soften length
//
// Return      :  External potential
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real ExternalPot( const double x, const double y, const double z, const double Time, const double UserArray[] )
{

   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
   const real   GM     = (real)UserArray[3];
   const real   dx     = (real)(x - Cen[0]);
   const real   dy     = (real)(y - Cen[1]);
   const real   dz     = (real)(z - Cen[2]);
   const real   _r     = 1.0/SQRT( dx*dx + dy*dy + dz*dz );

   return -GM*_r;

} // FUNCTION : ExternalPot



#endif // #ifdef GRAVITY



#endif // #ifndef __EXTERNALPOT__
