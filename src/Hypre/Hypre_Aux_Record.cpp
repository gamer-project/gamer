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
         fprintf( File, "# ====================================================================================================\n");
         fprintf( File, "# Hypre package information\n" );
         fprintf( File, "# ====================================================================================================\n");
#        ifdef HYPRE_RELEASE_DATE
         fprintf( File, "# %-30s : %s\n", "HYPRE_RELEASE_DATE",    HYPRE_RELEASE_DATE );
#        endif
#        ifdef HYPRE_RELEASE_NUMBER
         fprintf( File, "# %-30s : %d\n", "HYPRE_RELEASE_NUMBER",  HYPRE_RELEASE_NUMBER );
#        endif
#        ifdef HYPRE_RELEASE_VERSION
         fprintf( File, "# %-30s : %s\n", "HYPRE_RELEASE_VERSION", HYPRE_RELEASE_VERSION );
#        endif

#        ifdef HYPRE_DEBUG
         fprintf( File, "# %-30s : %d\n", "HYPRE_DEBUG",           HYPRE_DEBUG );
#        else
         fprintf( File, "# %-30s : %d\n", "HYPRE_DEBUG",           0 );
#        endif
#        ifdef HYPRE_HAVE_MPI
         fprintf( File, "# %-30s : %d\n", "HYPRE_HAVE_MPI",        HYPRE_HAVE_MPI );
#        else
         fprintf( File, "# %-30s : %d\n", "HYPRE_HAVE_MPI",        0 );
#        endif
#        ifdef HYPRE_USING_OPENMP
         fprintf( File, "# %-30s : %d\n", "HYPRE_USING_OPENMP",    HYPRE_USING_OPENMP );
#        else
         fprintf( File, "# %-30s : %d\n", "HYPRE_USING_OPENMP",    0 );
#        endif

         fprintf( File, "# ====================================================================================================\n");
         fprintf( File, "# Hypre runtime parameters\n" );
         fprintf( File, "# ====================================================================================================\n");
         fprintf( File, "# Maximum iteration:  %d\n",      HYPRE_MAX_ITER );
         fprintf( File, "# Relative tolerance: %22.16e\n", HYPRE_REL_TOL );
         fprintf( File, "# Absolute tolerance: %22.16e\n", HYPRE_ABS_TOL );
         fprintf( File, "# ====================================================================================================\n");
         fprintf( File, "#%19s  %5s  %10s  %14s", "Solver", "Lv", "Step", "Counter" );
         fprintf( File, "  %20s  %24s", "Iteration", "Residual" );
         fprintf( File, "\n" );
         fclose( File );
      }

      FirstTime = false;
   } // if ( FirstTime )

   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );
      fprintf( File, "%20s  %5d  %10ld  %14ld  %20d  %24.16e\n", SolverName, lv, Step, AdvanceCounter[lv], iteration, residual );
      fclose( File );
   }

} // FUNCTION : Hypre_Aux_Record



#endif // #ifdef SUPPORT_HYPRE
