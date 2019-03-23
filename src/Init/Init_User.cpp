#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_User();

// this function pointer may be overwritten by various test problem initializers
void (*Init_User_Ptr)() = Init_User;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User
// Description :  Add user-defined initialization
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_User()
{

} // FUNCTION : Init_User
