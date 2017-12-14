#include "GAMER.h"

#ifdef SUPPORT_LIBYT




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_Init
// Description :  Initialize the yt inline analysis
//
// Note        :  1. This function must be invoked once and only once during the entire simulation 
//
// Parameter   :  argc : Argument count
//                argv : Argument vector
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_Init( int argc, char *argv[] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// libyt runtime parameters
   yt_param_libyt param_libyt;

// verbose level
   param_libyt.verbose = YT_VERBOSE;

// YT analysis script without the ".py" extension (default="yt_inline_script")
   param_libyt.script  = YT_SCRIPT;

// initialize libyt
   if ( yt_init( argc, argv, &param_libyt ) != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_init() failed !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_Init



#endif // #ifdef SUPPORT_LIBYT
