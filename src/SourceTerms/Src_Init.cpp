#include "GAMER.h"



// prototypes of built-in source terms
void Src_Init_Deleptonization();

// this function pointer can be set by a test problem initializer for a non-built-in source term
void (*Src_Init_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Init
// Description :  Initialize the source terms
//
// Note        :  1. Invoked by Init_GAMER()
//                2. For a non-built-in source term, "Src_Init_User_Ptr" must be set by a test problem initializer
//                   in advance
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Src_Init()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check if at least one source term is activated
   if ( SRC_DELEPTONIZATION  ||  SRC_USER )  Src_Any = true;
   else                                      Src_Any = false;


// initialize source terms
   /*
   if ( SRC_DELEPTONIZATION )
      Src_Init_Deleptonization();
      */

   if ( SRC_USER )
   {
      if ( Src_Init_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "Src_Init_User_Ptr == NULL !!\n" );

      Src_Init_User_Ptr();
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Src_Init
