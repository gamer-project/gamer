#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  FindFather
// Description :  Construct the patch relation between levels "lv" and "lv-1"
//
// Note        :  a. This function assumes that the relations between levels "0, 1, 2, ..., lv-1" have already
//                   been constructed correctly
//                b. Currently this function only works for the functions "Init_ByRestart_v1, Init_ByRestart_v2,
//                   and Init_ByRestart_HDF5"
//                c. Only work on the "real" patches
//
// Parameter   :  lv    : Target refinement level
//                Mode  : 1 --> Find the father patch hierarchically from the base level
//                        2 --> Find the father patch by searching over all patches at level "lv-1"
//-------------------------------------------------------------------------------------------------------
void FindFather( const int lv, const int Mode )
{

   if ( lv <= 0  ||  lv > NLEVEL-1 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if ( Mode < 1  ||  Mode > 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Mode", Mode );


   const int scale0      = amr->scale[0];
   const int PS0         = PATCH_SIZE*scale0;
   const int BaseP_Disp  = 2;
   const int NBaseP[3]   = {  NX0[0]/PATCH_SIZE + 4,
                              NX0[1]/PATCH_SIZE + 4,
                              NX0[2]/PATCH_SIZE + 4  };
   const int RankDisp[3] = {  MPI_Rank_X[0]*NX0[0]*scale0,
                              MPI_Rank_X[1]*NX0[1]*scale0,
                              MPI_Rank_X[2]*NX0[2]*scale0  };

   int BaseP_ID, BaseP_xyz[3], GrandPaLv, GrandPaPID, FaPS, LocalPos[3], LocalID_1D;
   int FaPID = -1, LocalID = -1;
   int *Corner = NULL, *GrandPaCorner = NULL, *FaCorner = NULL;


   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
      Corner = amr->patch[0][lv][PID0]->corner;

      if ( Mode == 1 ) // construct the relation hierarchically from the base level
      {
//       a. find the ancestor patch at the base level
         for (int d=0; d<3; d++)    BaseP_xyz[d] = ( Corner[d] - RankDisp[d] ) / PS0  +  BaseP_Disp;

         BaseP_ID = BaseP_xyz[2]*NBaseP[1]*NBaseP[0] + BaseP_xyz[1]*NBaseP[0] + BaseP_xyz[0];
         FaPID    = BaseP[ BaseP_ID ];


//       b. find the father patch
         for (int FaLv=1; FaLv<lv; FaLv++)
         {
            GrandPaLv     = FaLv - 1;
            GrandPaPID    = FaPID;
            GrandPaCorner = amr->patch[0][GrandPaLv][GrandPaPID]->corner;
            FaPS          = PATCH_SIZE*amr->scale[FaLv];

            for (int d=0; d<3; d++)
            LocalPos[d]   = ( Corner[d] - GrandPaCorner[d] ) / FaPS;

            LocalID_1D    = LocalPos[2]*4 + LocalPos[1]*2 + LocalPos[0];

            switch ( LocalID_1D )
            {
               case 0:  LocalID = 0;   break;
               case 1:  LocalID = 1;   break;
               case 2:  LocalID = 2;   break;
               case 3:  LocalID = 4;   break;
               case 4:  LocalID = 3;   break;
               case 5:  LocalID = 6;   break;
               case 6:  LocalID = 5;   break;
               case 7:  LocalID = 7;   break;
               default:
                  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "LocalID_1D", LocalID_1D );
            }

            FaPID = amr->patch[0][GrandPaLv][GrandPaPID]->son + LocalID;

         } // for (int FaLv=1; FaLv<lv; FaLv++)
      } // if ( Mode == 1 )


      else // construct relation by searching over all patches at "lv-1"
      {
         for (int FaPID_Candidate=0; FaPID_Candidate<amr->NPatchComma[lv-1][1]; FaPID_Candidate++)
         {
            FaCorner = amr->patch[0][lv-1][FaPID_Candidate]->corner;

            if ( Corner[0] == FaCorner[0]  &&  Corner[1] == FaCorner[1]  &&  Corner[2] == FaCorner[2] )
            {
               FaPID = FaPID_Candidate;
               break;
            }

#           ifdef GAMER_DEBUG
            if ( FaPID_Candidate == amr->NPatchComma[lv-1][1]-1 )
               Aux_Error( ERROR_INFO, "No FaPID_Candidate is found !!\n" );
#           endif
         }
      } // if ( Mode == 1 ) ... else ...


//    store the relation
      amr->patch[0][lv-1][FaPID]->son = PID0;

      for (int PID=PID0; PID<PID0+8; PID++)  amr->patch[0][lv][PID]->father = FaPID;

   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

} // FUNCTION : FindFather
