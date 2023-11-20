#ifndef __CUFLU_MICROSHARED__
#define __CUFLU_MICROSHARED__


#include "CUFLU.h"

//-----------------------------------------------------------------------------------------
// Function    : MC_limiter
// Description : Monotonized central (MC) slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
static real MC_limiter( const real a, const real b )
{

   return minmod( (real)2.0*minmod(a,b), (real)0.5*(a+b) );

} // FUNCTION : MC_limiter


//-----------------------------------------------------------------------------------------
// Function    : minmod
// Description : Minmod slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
static real minmod( const real a, const real b )
{

   if      ( a > (real)0.0  &&  b > (real)0.0 )    return FMIN(a, b);
   else if ( a < (real)0.0  &&  b < (real)0.0 )    return FMAX(a, b);
   else                                            return (real)0.0;

} // FUNCTION : minmod

#endif // #ifndef __CUFLU_MICROSHARED__