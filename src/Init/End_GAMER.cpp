#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  End_GAMER
// Description :  Put everything you want to do before the end of program right here
//-------------------------------------------------------------------------------------------------------
void End_GAMER()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... \n", __FUNCTION__ );


#  ifdef TIMING
   Aux_DeleteTimer();
#  endif

   End_MemFree();

   End_TestProb();

#  ifdef GRAVITY
   End_FFTW();
#  endif

#  ifdef SUPPORT_LIBYT
   YT_End();
#  endif


   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
      Aux_Message( stdout, "\n\n~ GAME OVER ~\n\n\n" );
   }

   MPI_Finalize();

   exit( 0 );

} // FUNCTION : End_GAMER
