#include "GAMER.h"

extern void (*End_User_Ptr)();




//-------------------------------------------------------------------------------------------------------
// Function    :  End_GAMER
// Description :  Put everything you want to do before terminating the program right here
//
// Note        :  1. The function pointer "End_User_Ptr" points to "End_User()" by default
//                   but may be overwritten by various test problem initializers
//-------------------------------------------------------------------------------------------------------
void End_GAMER()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


#  ifdef TIMING
   Aux_DeleteTimer();
#  endif

   End_MemFree();

   if ( End_User_Ptr != NULL )   End_User_Ptr();

#  ifdef GRAVITY
   End_FFTW();
#  endif

#  ifdef SUPPORT_LIBYT
   YT_End();
#  endif


#  ifdef SUPPORT_GRACKLE
   Grackle_End();
#  endif


   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
      Aux_Message( stdout, "\n\n~ GAME OVER ~\n\n\n" );
   }

   MPI_Finalize();

   exit( 0 );

} // FUNCTION : End_GAMER
