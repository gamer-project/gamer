#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void End_User_Template();

// this function pointer must be set by a test problem initializer
void (*End_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  End_User_Template
// Description :  Template of user-specified operations before terminating the program
//
// Note        :  1. Invoked by End_GAMER() using the function pointer "End_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_User_Template()
{

} // FUNCTION : End_User_Template
