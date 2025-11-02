#ifndef __CUFLU_MICROSHARED__
#define __CUFLU_MICROSHARED__

#include "CUFLU.h"

#if defined( VISCOSITY ) || defined( CONDUCTION ) || defined( CR_DIFFUSION )

//-----------------------------------------------------------------------------------------
// Function    : minmod
// Description : Minmod slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real minmod( const real a, const real b )
{

   if      ( a > (real)0.0  &&  b > (real)0.0 )    return FMIN(a, b);
   else if ( a < (real)0.0  &&  b < (real)0.0 )    return FMAX(a, b);
   else                                            return (real)0.0;

} // FUNCTION : minmod

//-----------------------------------------------------------------------------------------
// Function    : MC_limiter
// Description : Monotonized central (MC) slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real MC_limiter( const real a, const real b )
{

   return minmod( (real)2.0*minmod(a,b), (real)0.5*(a+b) );

} // FUNCTION : MC_limiter

//-----------------------------------------------------------------------------------------
// Function    : van_leer 
// Description : vanLeer slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real van_leer( const real a, const real b )
{

   return ( a*b > (real)0.0 ) ? ( (real)2.0*a*b / (a+b) ) : (real)0.0;

} // FUNCTION : van_leer

//-----------------------------------------------------------------------------------------
// Function    : van_leer2 
// Description : vanLeer slope limiter called twice
//
// Parameter   : a, b, c, d : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real van_leer2( const real a, const real b, const real c, const real d )
{

   return van_leer( van_leer(a,b), van_leer(c,d) );

} // FUNCTION : van_leer2

#endif // #if defined( VISCOSITY ) || defined( CONDUCTION ) || defined( CR_DIFFUSION )

#endif // #ifndef __CUFLU_MICROSHARED__
