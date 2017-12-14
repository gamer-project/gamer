#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Output_User();

// this function pointer may be overwritten by various test problem initializers
void (*Output_User_Ptr)() = Output_User;




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_User
// Description :  User-specified output routine
//
// Note        :  1. Invoked by "Output_DumpData" using the function pointer "Output_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Enabled by the runtime option "OPT__OUTPUT_USER"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Output_User()
{

} // FUNCTION : Output_User
