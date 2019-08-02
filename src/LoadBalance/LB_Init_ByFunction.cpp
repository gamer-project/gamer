#include "GAMER.h"

#ifdef LOAD_BALANCE

void Init_ByFunction_AssignData( const int lv );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Init_ByFunction
// Description :  Set up the initial condition by invoking Init_ByFunction_AssignData()
//
// Note        :  1. Alternative function of Init_ByFunction()
//                2. Optimize load balancing on a level-by-level basis
//                   --> By invoking LB_Init_LoadBalance()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LB_Init_ByFunction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const bool   FindHomePatchForPar_Yes = true;
   const bool   Redistribute_Yes        = true;
   const bool   ResetLB_Yes             = true;
#  ifdef PARTICLE
   const double Par_Weight              = amr->LB->Par_Weight;
#  else
   const double Par_Weight              = 0.0;
#  endif


// construct all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing level %d ...\n", lv );

//    initialize the base level
      if ( lv == 0 )
         Init_UniformGrid( lv, FindHomePatchForPar_Yes );

//    initialize the refinement levels
      else
      {
         Flag_Real( lv-1, USELB_YES );

         LB_Init_Refine( lv-1 );
      }

//    get the total number of real patches
      Mis_GetTotalPatchNumber( lv );

//    assign initial condition on grids
      Init_ByFunction_AssignData( lv );

//    load balance
      LB_Init_LoadBalance( Redistribute_Yes, Par_Weight, ResetLB_Yes, lv );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing level %d ... done\n", lv );

   } // for (int lv=0; lv<NLEVEL; lv++)


// restrict all variables to be consistent with the finite volume scheme
   if ( OPT__INIT_RESTRICT )
   for (int lv=NLEVEL-2; lv>=0; lv--)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Restricting level %d ... ", lv );

      Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], amr->MagSg[lv+1], amr->MagSg[lv], NULL_INT, NULL_INT, _TOTAL, _MAG );

      LB_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_RESTRICT, _TOTAL, _MAG, NULL_INT );

      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_YES );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : LB_Init_ByFunction



#endif // #ifdef LOAD_BALANCE
