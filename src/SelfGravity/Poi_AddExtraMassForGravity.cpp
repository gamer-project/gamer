#include "GAMER.h"

#ifdef GRAVITY


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static real Poi_AddExtraMassForGravity( const double x, const double y, const double z, const double Time,
                                        const int lv, double AuxArray[] );

// this function pointer may be overwritten by various test problem initializers
real (*Poi_AddExtraMassForGravity_Ptr)( const double x, const double y, const double z, const double Time,
                                        const int lv, double AuxArray[] ) = Poi_AddExtraMassForGravity;




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_AddExtraMassForGravity
// Description :  Add extra mass source when computing gravity
//
// Note        :  1. Invoked by several functions using the function pointer "Poi_AddExtraMassForGravity_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Mass introduced here will only be used for computing gravity (i.e., when solving the Poisson eq.)
//                   --> It will NOT be used when solving other eqs. (e.g., hydro/MHD/Schroedinger)
//                   --> It will NOT be stored in the output data
//                   --> It WILL be included when computing the average density in Poi_GetAverageDensity(),
//                       which will then be used as DC in the periodic Poisson solver
//
// Parameter   :  x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  Mass density at (x, y, z, t)
//-------------------------------------------------------------------------------------------------------
real Poi_AddExtraMassForGravity( const double x, const double y, const double z, const double Time,
                                 const int lv, double AuxArray[] )
{

   real density = (real)0.0;
   return density;

} // FUNCTION : Poi_AddExtraMassForGravity



#endif // #ifdef GRAVITY
