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

// checks
#  if ( !defined SERIAL  &&  HYPRE_HAVE_MPI != 1 )
   Aux_Error( ERROR_INFO, "GAMER MPI (enabled) != HYPRE MPI (disable) !!\n" );
#  endif

#  if ( defined SERIAL  &&  HYPRE_HAVE_MPI == 1 )
   Aux_Error( ERROR_INFO, "GAMER MPI (disabled) != HYPRE MPI (enable) !!\n" );
#  endif

#  if ( defined GPU  &&  HYPRE_USING_CUDA != 1 )
   Aux_Error( ERROR_INFO, "GAMER GPU (enabled) != HYPRE GPU (disable) !!\n" );
#  endif

#  if ( !defined GPU  &&  HYPRE_USING_CUDA == 1 )
   Aux_Error( ERROR_INFO, "GAMER GPU (disabled) != HYPRE GPU (enable) !!\n" );
#  endif


// warning
#  if ( defined OPENMP  &&  HYPRE_USING_OPENMP != 1 )
   Aux_Message( stdout, "WARNING : GAMER OPENMP (enabled) != HYPRE OPENMP (disable) !!\n" );
#  endif

#  if ( !defined OPENMP  &&  HYPRE_USING_OPENMP == 1 )
   Aux_Message( stdout, "WARNING : GAMER OPENMP (disabled) != HYPRE OPENMP (enable) !!\n" );
#  endif


// initialize Hypre
   HYPRE_CHECK_FUNC(   HYPRE_Initialize()   );

#  ifdef GPU
   for (int r=0; r<MPI_NRank; r++)
   {
      if ( MPI_Rank == r )
      {
         printf( "GPU Info of Rank %d\n", MPI_Rank );
         HYPRE_PrintDeviceInfo();
      }
      MPI_Barrier( MPI_COMM_WORLD );
   }

   HYPRE_CHECK_FUNC(   HYPRE_SetExecutionPolicy( HYPRE_EXEC_DEVICE )   );
   HYPRE_CHECK_FUNC(   HYPRE_SetMemoryLocation( HYPRE_MEMORY_DEVICE )   );
#  endif

// TODO : print hypre info? check hypre info

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Hypre_Init
#endif // #ifdef SUPPORT_HYPRE
