#include "GAMER.h"



// prototypes of built-in source terms
#if ( MODEL == HYDRO )
void Src_End_Deleptonization();
void Src_End_ExactCooling();
#endif

// this function pointer can be set by a test problem initializer for a non-built-in source term
void (*Src_End_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_End
// Description :  Free the resources used by the source-term routines
//
// Note        :  1. Invoked by End_GAMER()
//                2. Set "Src_End_User_Ptr" in a test problem initializer for a non-built-in source term
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Src_End()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// stop all source terms
#  if ( MODEL == HYDRO )
   if ( SrcTerms.Deleptonization )
      Src_End_Deleptonization();
   if ( SrcTerms.ExactCooling )
      Src_End_ExactCooling(); 
#  endif

// users may not define Src_End_User_Ptr
   if ( SrcTerms.User  &&  Src_End_User_Ptr )
      Src_End_User_Ptr();


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Src_End
