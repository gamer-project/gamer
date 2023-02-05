#include "GAMER.h"


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Output_UserWorkBeforeOutput_Template();

// to enable this feature, set this function pointer by a test problem initializer
void (*Output_UserWorkBeforeOutput_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_UserWorkBeforeOutput_Template
// Description :  Template of user-specified work before dumping data
//
// Note        :  1. Invoked by Output_DumpData() using the function pointer "Output_UserWorkBeforeOutput_Ptr"
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Output_UserWorkBeforeOutput_Template()
{


} // FUNCTION : Output_UserWorkBeforeOutput_Template
