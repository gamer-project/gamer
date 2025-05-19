#include "GAMER.h"




#ifdef SUPPORT_HYPRE
//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_Init
// Description :  Initialize Hypre
//
// Note        :  1. Invoked by Init_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Hypre_Init()
{

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// initialize Hypre
   HYPRE_CHECK_FUNC(   HYPRE_Initialize()   );

// TODO : print hypre info? check hypre info

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Hypre_Init
#endif // #ifdef SUPPORT_HYPRE
