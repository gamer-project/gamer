#include "GAMER.h"




#ifdef SUPPORT_HYPRE
static void Hypre_Solve_SStructSysPFMG( int *N_iter, real *final_res_norm );
static void Hypre_Solve_SStructSplit  ( int *N_iter, real *final_res_norm );



void Hypre_Solve( const Hypre_Solver_t Solver, int *N_iter, real *final_res_norm )
{

   switch ( Solver )
   {
      case HYPRE_SOLVER_SSTRUCT_SYS_PFMG: Hypre_Solve_SStructSysPFMG( N_iter, final_res_norm );  break;
      case HYPRE_SOLVER_SSTRUCT_SPLIT:    Hypre_Solve_SStructSplit  ( N_iter, final_res_norm );  break;
      default: Aux_Error( ERROR_INFO, "Unknown HYPRE_SOLVER: %d !!\n", HYPRE_SOLVER );
   } // switch ( Solver )

#  ifdef DEBUG_HYPRE
   if ( MPI_Rank == 0 )   Aux_Message( stdout, "%s: Iteration: %d, Residual norm: %24.16e\n", __FUNCTION__, N_iter, final_res_norm );
#  endif
   // TODO : warning message for reaching the max iter

} // FUNCTION : Hypre_Solve



void Hypre_Solve_SStructSysPFMG( int *N_iter, real *final_res_norm )
{

   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGCreate( HYPRE_MPI_COMM, &Hypre_solver )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGSetPrintLevel( Hypre_solver, HYPRE_PRINT_LEVEL )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGSetLogging( Hypre_solver, HYPRE_ENABLE_LOGGING )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGSetZeroGuess( Hypre_solver )   );

// Set sysPFMG parameters
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGSetTol( Hypre_solver, HYPRE_REL_TOL )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGSetMaxIter( Hypre_solver, HYPRE_MAX_ITER )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGSetNumPreRelax( Hypre_solver, HYPRE_NPRE_RELAX )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGSetNumPostRelax( Hypre_solver, HYPRE_NPOST_RELAX )   );

// do the setup
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGSetup( Hypre_solver, Hypre_A, Hypre_b, Hypre_x )   );

// do the solve
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGSolve( Hypre_solver, Hypre_A, Hypre_b, Hypre_x )   );

// get some info
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm( Hypre_solver, final_res_norm )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGGetNumIterations( Hypre_solver, N_iter )   );

// clean up
   HYPRE_CHECK_FUNC(   HYPRE_SStructSysPFMGDestroy( Hypre_solver )   );

} // FUNCTION : Hypre_Solve_SStructSysPFMG



void Hypre_Solve_SStructSplit( int *N_iter, real *final_res_norm )
{

   HYPRE_CHECK_FUNC(   HYPRE_SStructSplitCreate( HYPRE_MPI_COMM, &Hypre_solver )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructSplitSetZeroGuess( Hypre_solver )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructSplitSetMaxIter( Hypre_solver, HYPRE_MAX_ITER )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructSplitSetTol( Hypre_solver, HYPRE_REL_TOL )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructSplitSetup( Hypre_solver, Hypre_A, Hypre_b, Hypre_x )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructSplitSolve( Hypre_solver, Hypre_A, Hypre_b, Hypre_x )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructSplitGetFinalRelativeResidualNorm( Hypre_solver, final_res_norm )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructSplitGetNumIterations( Hypre_solver, N_iter )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructSplitDestroy( Hypre_solver )   );

} // FUNCTION : Hypre_Solve_SStructSplit
#endif // #ifdef SUPPORT_HYPRE
