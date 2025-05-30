#include "GAMER.h"




#ifdef SUPPORT_HYPRE



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_Center
// Description :  Record various center coordinates
//
// Note        :  1.
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Hypre_Aux_Record( const char SolverName[], const int lv, const int iteration, const real residual )
{

   static bool FirstTime = true;
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__Hypre", OUTPUT_DIR );

   if ( FirstTime )
   {
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist( FileName ) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File = fopen( FileName, "a" );
         fprintf( File, "# Maximum iteration: %d\n", HYPRE_MAX_ITER );
         fprintf( File, "# Relative tolerance: %24.16e\n", HYPRE_REL_TOL );
         fprintf( File, "# Absolute tolerance: %24.16e\n", HYPRE_ABS_TOL );
         fprintf( File, "# ====================================================================================================\n");
         fprintf( File, "#%19s  %5s  %10s  %14s", "Solver", "Lv", "Step", "Counter" );
         fprintf( File, "  %20s  %24s", "Iteration", "Residual" );
         fprintf( File, "\n" );
         fclose( File );
      }

      FirstTime = false;
   }

   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );
      fprintf( File, "%20s  %5d  %10ld  %14ld  %20d  %+24.16e\n", SolverName, lv, Step, AdvanceCounter[lv], iteration, residual );
      fclose( File );
   }

} // FUNCTION : Hypre_Aux_Record



#endif // #ifdef SUPPORT_HYPRE
