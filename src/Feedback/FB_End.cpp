#include "GAMER.h"

#ifdef FEEDBACK



// prototypes of built-in feedbacks
void FB_End_SNE();


// user-specified feedback to be set by a test problem initializer
void (*FB_End_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End
// Description :  Free the resources used by the feedback routines
//
// Note        :  1. Invoked by End_GAMER()
//                2. Set "FB_End_User_Ptr" in a test problem initializer for a non-built-in feedback
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   if ( FB_USER  &&  FB_End_User_Ptr != NULL )  FB_End_User_Ptr();


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : FB_End



#endif // #ifdef FEEDBACK
