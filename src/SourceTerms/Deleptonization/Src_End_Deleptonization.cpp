#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_End_Deleptonization
// Description :  Free the resources used by the deleptonization source term
//
// Note        :  1. Invoked by Src_End()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Src_End_Deleptonization()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   %s ...\n", __FUNCTION__ );



   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   %s ... done\n", __FUNCTION__ );

} // FUNCTION : Src_End_Deleptonization
