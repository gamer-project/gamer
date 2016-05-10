#include "Copyright.h"
#include "GAMER.h"


void (*Output_TestProbErr_Ptr)( const bool BaseOnly ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_TestProbErr
// Description :  Compare and output the numerical and analytical solutions for the chosen test problem 
//
// Note        :  1. This function will invoke the function pointer "Output_TestProbErr_Ptr", which should be 
//                   specified in the function "Init_TestProb" of the target test problem
//                   (e.g., GAMER/test_problem/Model_ELBDM/Jeans_Instability/Init_TestProb.cpp). 
//                2. Not all test problems support this function
//
// Parameter   :  BaseOnly :  Only output the base-level data
//-------------------------------------------------------------------------------------------------------
void Output_TestProbErr( const bool BaseOnly )
{  

// output the test problem error if the target function is specified
   if ( Output_TestProbErr_Ptr != NULL )  Output_TestProbErr_Ptr( BaseOnly );

} // FUNCTION : Output_TestProbErr
