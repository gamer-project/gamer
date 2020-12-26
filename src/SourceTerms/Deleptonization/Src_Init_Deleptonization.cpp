#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Init_Deleptonization
// Description :  Initialize the deleptonization source term
//
// Note        :  1. Invoked by Src_Init()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Src_Init_Deleptonization()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   %s ...\n", __FUNCTION__ );



   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   %s ... done\n", __FUNCTION__ );

} // FUNCTION : Src_Init_Deleptonization
