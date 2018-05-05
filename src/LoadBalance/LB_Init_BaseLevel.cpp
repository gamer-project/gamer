#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Init_BaseLevel
// Description :  Construct the base level during initialization for LOAD_BALANCE
//
// Note        :  1. Invoked by LB_Init_ByFunction()
//                2. Invoke Par_FindHomePatch_Base() to associate particles with base-level patches
//                3. Invoke LB_Init_LoadBalance() to initialize load-balancing on the base level
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LB_Init_BaseLevel()
{

// 1. set up the load-balance cut points amr->LB->CutPoint[] on the base level
   const bool   InputLBIdx0AndLoad_Yes = true;
   const double ParWeight_Zero         = 0.0;
   const int    NPG_Total              = (NX0_TOT[0]/PS2)*(NX0_TOT[1]/PS2)*(NX0_TOT[2]/PS2);
   const int    scale0                 = amr->scale[0];
   const int    lv0                    = 0;

   long   *LBIdx0_AllRank = NULL;
   double *Load_AllRank   = NULL;
   int     Cr[3], counter;

// 1.1 prepare LBIdx and load-balance weighting of each **patch group** for LB_SetCutPoint()
   if ( MPI_Rank == 0 )
   {
      LBIdx0_AllRank = new long   [NPG_Total];
      Load_AllRank   = new double [NPG_Total];

      counter = 0;

      for (int k=0; k<NX0_TOT[2]/PS2; k++)   {  Cr[2] = k*PS2*scale0;
      for (int j=0; j<NX0_TOT[1]/PS2; j++)   {  Cr[1] = j*PS2*scale0;
      for (int i=0; i<NX0_TOT[0]/PS2; i++)   {  Cr[0] = i*PS2*scale0;

         LBIdx0_AllRank[counter]  = LB_Corner2Index( lv0, Cr, CHECK_ON );
         LBIdx0_AllRank[counter] -= LBIdx0_AllRank[counter] % 8;  // get the minimum LBIdx in each patch group
         Load_AllRank  [counter]  = 8.0;                          // assuming all patches have the same weighting == 1.0

         counter ++;
      }}}
   }

// 1.2 set CutPoint[]
//     --> do NOT consider load-balance weighting of particles since we have not assoicated particles with patches yet
   LB_SetCutPoint( lv0, NPG_Total, amr->LB->CutPoint[lv0], InputLBIdx0AndLoad_Yes, LBIdx0_AllRank, Load_AllRank,
                   ParWeight_Zero );

// 1.3 free memory
   if ( MPI_Rank == 0 )
   {
      delete [] LBIdx0_AllRank;
      delete [] Load_AllRank;
   }


// 2. allocate the real base-level patches
   const int PScale0 = PS1*scale0;

   for (int k=0; k<NX0_TOT[2]/PS2; k++)   {  Cr[2] = k*PS2*scale0;
   for (int j=0; j<NX0_TOT[1]/PS2; j++)   {  Cr[1] = j*PS2*scale0;
   for (int i=0; i<NX0_TOT[0]/PS2; i++)   {  Cr[0] = i*PS2*scale0;

      const long LBIdx0 = LB_Corner2Index( lv0, Cr, CHECK_ON );

      if (  LB_Index2Rank( lv0, LBIdx0, CHECK_ON ) == MPI_Rank  )
      {
         amr->pnew( lv0, Cr[0],         Cr[1],         Cr[2],         -1, true, true );
         amr->pnew( lv0, Cr[0]+PScale0, Cr[1],         Cr[2],         -1, true, true );
         amr->pnew( lv0, Cr[0],         Cr[1]+PScale0, Cr[2],         -1, true, true );
         amr->pnew( lv0, Cr[0],         Cr[1],         Cr[2]+PScale0, -1, true, true );
         amr->pnew( lv0, Cr[0]+PScale0, Cr[1]+PScale0, Cr[2],         -1, true, true );
         amr->pnew( lv0, Cr[0],         Cr[1]+PScale0, Cr[2]+PScale0, -1, true, true );
         amr->pnew( lv0, Cr[0]+PScale0, Cr[1],         Cr[2]+PScale0, -1, true, true );
         amr->pnew( lv0, Cr[0]+PScale0, Cr[1]+PScale0, Cr[2]+PScale0, -1, true, true );
      }

   }}}

   for (int m=1; m<28; m++)   amr->NPatchComma[lv0][m] = amr->num[lv0];


// 3. find the base-level home patches of all particles
#  ifdef PARTICLE
   Par_FindHomePatch_Base();
#  endif

} // FUNCTION : LB_Init_BaseLevel



#endif // #ifdef LOAD_BALANCE
