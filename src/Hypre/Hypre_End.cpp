#include "GAMER.h"




#ifdef SUPPORT_HYPRE
//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_End
// Description :  Free the resources used by the Hypre
//
// Note        :  1. Invoked by End_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Hypre_End()
{

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// Hypre_Free();
// finalize Hypre
   HYPRE_CHECK_FUNC(   HYPRE_Finalize()   );

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Hypre_End
#endif // #ifdef SUPPORT_HYPRE
