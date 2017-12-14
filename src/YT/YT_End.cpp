#include "GAMER.h"

#ifdef SUPPORT_LIBYT




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_End
// Description :  Finalize the yt inline analysis
//
// Note        :  1. This function must be invoked once and only once during the entire simulation
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_End()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   if ( yt_finalize() != YT_SUCCESS )   Aux_Error( ERROR_INFO, "yt_finalize() failed !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_End



#endif // #ifdef SUPPORT_LIBYT
