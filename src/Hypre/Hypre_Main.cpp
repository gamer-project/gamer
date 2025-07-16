#include "GAMER.h"



#ifdef SUPPORT_HYPRE
//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_Solver
// Description :
//
// Parameter   :
//-------------------------------------------------------------------------------------------------------
void Hypre_Solver( const Hypre_SolveType_t SolveType, const int lv, const double TimeNew, const double TimeOld,
                   const double dt_in, const double Poi_Coeff, const int SaveSg_Flu, const int SaveSg_Mag,
                   const int SaveSg_Pot )
{

// check
// TODO

   if ( NPatchTotal[lv] == 0 )   return;

   int  N_iter;
   real final_residual;

// 1. prepare single level
   Hypre_PrepareSingleLevel( SolveType, lv );

// 2. fill arrays
   Hypre_FillArrays( SolveType, lv, TimeNew, Poi_Coeff, SaveSg_Pot );

// 3. solve
   Hypre_Solve( HYPRE_SOLVER, &N_iter, &final_residual );

// 4. collect the answer from all ranks
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGather( Hypre_x )   );

// 5. update gamer array
   Hypre_UpdateArrays( SolveType, lv, SaveSg_Flu, SaveSg_Mag, SaveSg_Pot );

// 6. free
   Hypre_Free();

// 7. record Hypre informaiton
   Hypre_Aux_Record( SolveType, lv, N_iter, final_residual );

} // FUNCTION : Hypre_Solver
#endif // #ifdef SUPPORT_HYPRE
