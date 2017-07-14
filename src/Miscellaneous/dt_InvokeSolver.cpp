#include "Copyright.h"
#include "GAMER.h"

double dt_min_for_solver;




//-------------------------------------------------------------------------------------------------------
// Function    :  dt_InvokeSolver
// Description :  Invoke the time-step solver
//
// Note        :  1. Invoked by Mis_GetTimeStep()
//                2. The global variable "dt_min_for_solver" will be set by dt_Close()
//
// Parameter   :  TSolver : Target dt solver
//                          --> DT_SOLVER_HYDRO_CFL     : HYDRO CFL
//                              DT_SOLVER_HYDRO_GRAVITY : HYDRO gravity
//                              DT_SOLVER_ELBDM_GRAVITY : ELBDM gravity
//                lv      : Target refinement level
//-------------------------------------------------------------------------------------------------------
double dt_InvokeSolver( const Solver_t TSolver, const int lv )
{

// initialize it as an extremely large value, which will be reset by dt_Close()
   dt_min_for_solver = HUGE_NUMBER;


// invoke the target dt solver
   InvokeSolver( TSolver, lv, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_INT, NULL_INT, false, false );


// get the minimum dt among all ranks
   double dt_min_all_rank;
   MPI_Allreduce( &dt_min_for_solver, &dt_min_all_rank, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );


   return dt_min_all_rank;

} // FUNCTION : dt_InvokeSolver
