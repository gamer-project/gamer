#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Microphysics_End
// Description :  Free the resources used by the microphysics routines
//
// Note        :  1. Invoked by End_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Microphysics_End()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// empty now

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Microphysics_End
