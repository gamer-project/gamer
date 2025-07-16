#include "GAMER.h"



#ifdef SUPPORT_HYPRE
static void Hypre_UpdateArrays_Poisson( const int lv, const int SaveSg_Pot );



//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_UpdateArrays
// Description :  TBF
//
// Note        :  TBF
//
// Parameter   :  TBF
//-------------------------------------------------------------------------------------------------------
void Hypre_UpdateArrays( const Hypre_SolveType_t SolveType, const int lv, const int SaveSg_Flu,
                         const int SaveSg_Mag, const int SaveSg_Pot )
{

   switch ( SolveType )
   {
#     ifdef GRAVITY
      case HYPRE_SOLVE_TYPE_POISSON:
         Hypre_UpdateArrays_Poisson( lv, SaveSg_Pot );
         break;
#     endif
      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SolveType", SolveType );
   } // switch ( SolveType )

} // FUNCTION : Hypre_UpdateArrays



#ifdef GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_UpdateArrays_Poisson
// Description :  TBF
//
// Note        :  TBF
//
// Parameter   :  TBF
//-------------------------------------------------------------------------------------------------------
void Hypre_UpdateArrays_Poisson( const int lv, const int SaveSg_Pot )
{

   const int part = 0; // single level only, no need to iterate parts
   const int var  = 0; // single variable only, no need to iterate variables

   real_hypre *pote;
#  ifdef GPU
   cudaMallocManaged( &pote, sizeof(real_hypre) * CUBE(PS1), cudaMemAttachGlobal );
#  else
   pote = new real_hypre [CUBE(PS1)];
#  endif
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         const int idx = IDX321( i, j, k, PS1, PS1 );
         amr->patch[ SaveSg_Pot ][lv][PID]->pot[k][j][i] = (real)pote[idx];
      } // i, j, k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

#  ifdef GPU
   cudaFree( pote );
#  else
   delete [] pote;
#  endif

} // FUNCTION : Hypre_UpdateArrays_Poisson
#endif // #ifdef GRAVITY



#endif // #ifdef SUPPORT_HYPRE
