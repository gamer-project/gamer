#include "GAMER.h"

// declared in dt_InvokeSolver.cpp
extern double dt_min_for_solver;




//-------------------------------------------------------------------------------------------------------
// Function    :  dt_Close
// Description :  Closing step for the dt solver
//
// Note        :  1. Get the minimum dt in all target patches
//                2. Store the minimum dt in the global variable "dt_min_for_solver" declared in dt_InvokeSolver.cpp
//                   --> "dt_min_for_solver" must be initialized as an extremely large value in advance
//
// Parameter   :  h_dt_Array_T : Host array storing the minimum dt in each target patch
//                NPG          : Number of target patch groups
//-------------------------------------------------------------------------------------------------------
void dt_Close( const real h_dt_Array_T[], const int NPG )
{

   const int NP = 8*NPG;

   for (int t=0; t<NP; t++)   dt_min_for_solver = fmin( dt_min_for_solver, (double)h_dt_Array_T[t] );

} // FUNCTION : dt_Close
