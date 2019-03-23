#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_UniformGrid
// Description :  Construct a uniform grid on a target level during initialization
//
// Note        :  1. After invoking this function, the target level will be fully occupied by patches
//                   --> As if "lv-1" is fully refined
//                2. Invoked by LB_Init_ByFunction() and Init_ByFile()
//                3. If FindHomePatchForPar == true --> invoke Par_FindHomePatch_UniformGrid() to associate
//                   particles with home patches on the target level
//                4. Work for both SERIAL and LOAD_BALANCE
//                   --> For LOAD_BALANCE, this function will invoke LB_SetCutPoint() to set the load-balance
//                       cut points and distribute patches accordingly
//                       --> But particle weighting is not taken into account
//                       --> One can invoke LB_Init_LoadBalance() later to optimize load-balancing with
//                           particle weighting
//                5. Do not allocate any buffer patch
//
// Parameter   :  lv                  : Target level
//                FindHomePatchForPar : Find home patches on lv for all particles in the particle repository
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_UniformGrid( const int lv, const bool FindHomePatchForPar )
{

   const int NPG_EachDim[3] = { (NX0_TOT[0]/PS2)*(1<<lv), (NX0_TOT[1]/PS2)*(1<<lv), (NX0_TOT[2]/PS2)*(1<<lv) };
   const int scale          = amr->scale[lv];

   int Cr[3];


// 1. set up the load-balance cut points amr->LB->CutPoint[] on lv
#  ifdef LOAD_BALANCE
   const bool   InputLBIdx0AndLoad_Yes = true;
   const double ParWeight_Zero         = 0.0;
   const long   NPG_Total              = (long)NPG_EachDim[0]*(long)NPG_EachDim[1]*(long)NPG_EachDim[2];

   long   *LBIdx0_AllRank = NULL;
   double *Load_AllRank   = NULL;
   int     counter;

// 1.1 prepare LBIdx and load-balance weighting of each **patch group** for LB_SetCutPoint()
   if ( MPI_Rank == 0 )
   {
      LBIdx0_AllRank = new long   [NPG_Total];
      Load_AllRank   = new double [NPG_Total];

      counter = 0;

      for (int k=0; k<NPG_EachDim[2]; k++)   {  Cr[2] = k*PS2*scale;
      for (int j=0; j<NPG_EachDim[1]; j++)   {  Cr[1] = j*PS2*scale;
      for (int i=0; i<NPG_EachDim[0]; i++)   {  Cr[0] = i*PS2*scale;

         LBIdx0_AllRank[counter]  = LB_Corner2Index( lv, Cr, CHECK_ON );
         LBIdx0_AllRank[counter] -= LBIdx0_AllRank[counter] % 8;  // get the minimum LBIdx in each patch group
         Load_AllRank  [counter]  = 8.0;                          // assuming all patches have the same weighting == 1.0

         counter ++;
      }}}
   }

// 1.2 set CutPoint[]
//     --> do NOT consider load-balance weighting of particles since we have not assoicated particles with patches yet
   LB_SetCutPoint( lv, NPG_Total, amr->LB->CutPoint[lv], InputLBIdx0AndLoad_Yes, LBIdx0_AllRank, Load_AllRank,
                   ParWeight_Zero );

// 1.3 free memory
   if ( MPI_Rank == 0 )
   {
      delete [] LBIdx0_AllRank;
      delete [] Load_AllRank;
   }
#  endif // #ifdef LOAD_BALANCE


// 2. allocate real patches on lv
   const int PScale = PS1*scale;

   for (int k=0; k<NPG_EachDim[2]; k++)   {  Cr[2] = k*PS2*scale;
   for (int j=0; j<NPG_EachDim[1]; j++)   {  Cr[1] = j*PS2*scale;
   for (int i=0; i<NPG_EachDim[0]; i++)   {  Cr[0] = i*PS2*scale;

#     ifdef LOAD_BALANCE
      const long LBIdx0 = LB_Corner2Index( lv, Cr, CHECK_ON );
      if (  LB_Index2Rank( lv, LBIdx0, CHECK_ON ) == MPI_Rank  )
#     endif
      {
         amr->pnew( lv, Cr[0],        Cr[1],        Cr[2],        -1, true, true );
         amr->pnew( lv, Cr[0]+PScale, Cr[1],        Cr[2],        -1, true, true );
         amr->pnew( lv, Cr[0],        Cr[1]+PScale, Cr[2],        -1, true, true );
         amr->pnew( lv, Cr[0],        Cr[1],        Cr[2]+PScale, -1, true, true );
         amr->pnew( lv, Cr[0]+PScale, Cr[1]+PScale, Cr[2],        -1, true, true );
         amr->pnew( lv, Cr[0],        Cr[1]+PScale, Cr[2]+PScale, -1, true, true );
         amr->pnew( lv, Cr[0]+PScale, Cr[1],        Cr[2]+PScale, -1, true, true );
         amr->pnew( lv, Cr[0]+PScale, Cr[1]+PScale, Cr[2]+PScale, -1, true, true );
      }

   }}}

   for (int m=1; m<28; m++)   amr->NPatchComma[lv][m] = amr->num[lv];


// 3. find the home patches on lv for all particles
#  ifdef PARTICLE
   if ( FindHomePatchForPar )
   {
      const bool OldParOnly_Yes = true;
      Par_FindHomePatch_UniformGrid( lv, OldParOnly_Yes, NULL_INT, NULL );
   }
#  endif

} // FUNCTION : Init_UniformGrid
