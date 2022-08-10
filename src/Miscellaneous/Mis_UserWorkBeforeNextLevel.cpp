#include "GAMER.h"


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Mis_UserWorkBeforeNextLevel_Template( const int lv, const double TimeNew, const double TimeOld, const double dt );

// to enable this feature, set this function pointer by a test problem initializer
void (*Mis_UserWorkBeforeNextLevel_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_UserWorkBeforeNextLevel_Template
// Description :  Template of user-specified work before entering the next AMR level in EvolveLevel()
//
// Note        :  1. Invoked by EvolveLevel() using the function pointer "Mis_UserWorkBeforeNextLevel_Ptr"
//
// Parameter   :  lv      : Target refinement level
//                TimeNew : Target physical time to reach
//                TimeOld : Physical time before update
//                dt      : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//-------------------------------------------------------------------------------------------------------
void Mis_UserWorkBeforeNextLevel_Template( const int lv, const double TimeNew, const double TimeOld, const double dt )
{


} // FUNCTION : Mis_UserWorkBeforeNextLevel_Template
