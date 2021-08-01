#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_User_Template();

// this function pointer must be set by a test problem initializer
void (*Init_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User_Template
// Description :  Template of user-defined initialization
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_User_Template()
{

} // FUNCTION : Init_User_Template
