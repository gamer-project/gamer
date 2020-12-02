#include "GAMER.h"

#ifdef GRAVITY


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Poi_UserWorkBeforePoisson_Template( const double Time, const int lv );

// to enable this feature, set this function pointer by a test problem initializer
void (*Poi_UserWorkBeforePoisson_Ptr)( const double Time, const int lv ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_UserWorkBeforePoisson_Template
// Description :  Template of user-specified work before invoking the Poisson solver
//
// Note        :  1. Invoked by Gra_AdvanceDt() using the function pointer "Poi_UserWorkBeforePoisson_Ptr"
//
// Parameter   :  Time : Target physical time
//                lv   : Target refinement level
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Poi_UserWorkBeforePoisson_Template( const double Time, const int lv )
{

   /*
// example: reset the auxiliary arrays of external acceleration and potential
//          --> useful when they are functions of time
   if ( OPT__EXT_ACC )
   SetExtAccAuxArray_PointMass( ExtAcc_AuxArray );

   if ( OPT__EXT_POT )
   SetExtPotAuxArray_PointMass( ExtPot_AuxArray );

#  ifdef GPU
   CUAPI_SetConstMemory_ExtAccPot();
#  endif
   */

} // FUNCTION : Poi_UserWorkBeforePoisson_Template



#endif // #ifdef GRAVITY
