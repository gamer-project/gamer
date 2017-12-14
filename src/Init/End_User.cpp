#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void End_User();

// this function pointer may be overwritten by various test problem initializers
void (*End_User_Ptr)() = End_User;




//-------------------------------------------------------------------------------------------------------
// Function    :  End_User
// Description :  User-specified operations before terminating the program
//
// Note        :  1. Invoked by "End_GAMER" using the function pointer "End_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_User()
{

} // FUNCTION : End_User
