#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Output_User_Template();

// this function pointer must be set by a test problem initializer
void (*Output_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_User_Template
// Description :  User-specified output template
//
// Note        :  1. Invoked by Output_DumpData() using the function pointer "Output_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__OUTPUT_USER"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Output_User_Template()
{

} // FUNCTION : Output_User_Template
