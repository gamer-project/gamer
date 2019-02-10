#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_IsFinite
// Description :  Check whether the input floating-point value is finite
//
// Note        :  1. Definition of "finite" --> not NaN, Inf, -Inf
//                2. Alternative to the built-in function isfinite() (or std::isfinite()) which is less portable
//
// Parameter   :  x : Floating-point value to be checked
//
// Return      :  1 : finite
//                0 : not finite
//-------------------------------------------------------------------------------------------------------
int Aux_IsFinite( const float x )
{

   if ( x != x  ||  x < -__FLT_MAX__  ||  x > __FLT_MAX__ )    return 0;
   else                                                        return 1;

} // FUNCTION : Aux_IsFinite



//-------------------------------------------------------------------------------------------------------
// Function overloading: double precision
//-------------------------------------------------------------------------------------------------------
int Aux_IsFinite( const double x )
{

   if ( x != x  ||  x < -__DBL_MAX__  ||  x > __DBL_MAX__ )    return 0;
   else                                                        return 1;

} // FUNCTION : Aux_IsFinite
