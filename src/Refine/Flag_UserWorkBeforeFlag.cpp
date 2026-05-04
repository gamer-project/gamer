#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Flag_UserWorkBeforeFlag_Template( const double Time, const int lv );

// to enable this feature, set this function pointer by a test problem initializer
void (*Flag_UserWorkBeforeFlag_Ptr)( const double Time, const int lv ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_UserWorkBeforeFlag_Template
// Description :  Template of user-specified work before flagging cells for refinement
//
// Note        :  1. Invoked by Flag_Real() using the function pointer "Flag_UserWorkBeforeFlag_Ptr"
//
// Parameter   :  Time : Target physical time
//                lv   : Target refinement level
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Flag_UserWorkBeforeFlag_Template( const double Time, const int lv )
{

} // FUNCTION : Flag_UserWorkBeforeFlag_Template
