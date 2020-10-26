#include "GAMER.h"

#ifdef GRAVITY


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static real Poi_AddExtraMassForGravity_Template( const double x, const double y, const double z, const double Time,
                                                 const int lv, double AuxArray[] );

// this function pointer must be set by a test problem initializer
real (*Poi_AddExtraMassForGravity_Ptr)( const double x, const double y, const double z, const double Time,
                                        const int lv, double AuxArray[] ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_AddExtraMassForGravity_Template
// Description :  Template of adding extra mass source when computing gravity
//
// Note        :  1. Invoked by several functions using the function pointer "Poi_AddExtraMassForGravity_Ptr",
//                   which must be set by a test problem initializer
//                2. Mass introduced here will only be used for computing gravity (i.e., when solving the Poisson eq.)
//                   --> It will NOT be used when solving other eqs. (e.g., hydro/MHD/Schroedinger)
//                   --> It will NOT be stored in the output data
//                   --> It WILL be included when computing the average density in Poi_GetAverageDensity(),
//                       which will then be used as DC in the periodic Poisson solver
//                3. Enabled by the runtime option "OPT__GRAVITY_EXTRA_MASS"
//
// Parameter   :  x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  Mass density at (x, y, z, t)
//-------------------------------------------------------------------------------------------------------
real Poi_AddExtraMassForGravity_Template( const double x, const double y, const double z, const double Time,
                                          const int lv, double AuxArray[] )
{

   real density = (real)0.0;
   return density;

} // FUNCTION : Poi_AddExtraMassForGravity_Template



#endif // #ifdef GRAVITY
