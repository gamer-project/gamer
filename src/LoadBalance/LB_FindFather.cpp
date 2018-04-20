#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_FindFather
// Description :  Construct the patch relation : son <-> father
//
// Note        :  1. LB_PaddedCr1DList and LB_PaddedCr1DList_IdxTable must be properly prepared
//                   at SonLv and FaLv
//                2. Father-buffer patches should be allocated in advance by "LB_AllocateBufferPatch_Father"
//                3. One should find father patches only for the "real" patches at SonLv (applying to
//                   sibling-buffer and father-buffer patches is not necessary)
//                   --> Father patches found by this function will never be "external" patches
//                4. Sibling-buffer and father-buffer patches at SonLv will have father indices == -1
//                5. For father patches with sons not home, their "father->son" relation should be set
//                   by "LB_FindSonNotHome"
//                6. This function does NOT initialize son indices as -1 for all patches at FaLv
//                7. This function does NOT initialize father indices as -1 for all patches at SonLv
//
// Parameter   :  SonLv          : Target refinement level of sons
//                SearchAllSon   : Search over all son patches at SonLv
//                NInput         : Number of target son patches (with LocalID==0) in "TargetSonPID0"
//                                 (useful only if SearchAllSon == false)
//                TargetSonPID0  : Lists recording all target son patches (with LocalID==0)
//                                 (useful only if SearchAllSon == false)
//-------------------------------------------------------------------------------------------------------
void LB_FindFather( const int SonLv, const bool SearchAllSon, const int NInput, int* TargetSonPID0 )
{

// nothing to do for the base level
   if ( SonLv == 0 )    return;


// check
   if ( SonLv < 0  ||  SonLv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonLv", SonLv );

   if ( !SearchAllSon  &&  NInput != 0  &&  TargetSonPID0 == NULL )
      Aux_Error( ERROR_INFO, "SonLv %d, NInput %d, TargetSonPID0 == NULL !!\n", SonLv, NInput );


   const int NTargetSon0 = ( SearchAllSon ) ? amr->NPatchComma[SonLv][1]/8 : NInput;
   const int FaLv        = SonLv - 1;

   ulong *Cr1D_Son0          = new ulong [NTargetSon0];
   int   *Cr1D_Son0_IdxTable = new int   [NTargetSon0];
   int   *Match_Son0         = new int   [NTargetSon0];

   int SonPID, SonPID0, FaID, FaPID;


// 1. nothing to do if there is no target real patch at SonLv
   if ( NTargetSon0 == 0 )
   {
      delete [] Cr1D_Son0;
      delete [] Cr1D_Son0_IdxTable;
      delete [] Match_Son0;

      return;
   }


// 2. construct the target son patch list
   if ( SearchAllSon )
   {
      TargetSonPID0 = new int [NTargetSon0];

      for (int t=0; t<NTargetSon0; t++)   TargetSonPID0[t] = t*8;
   }

#  ifdef GAMER_DEBUG
   for (int t=0; t<NTargetSon0; t++)
   {
      SonPID0 = TargetSonPID0[t];

      if ( SonPID0 >= amr->NPatchComma[SonLv][1] )
         Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 (%d) is not a real patch (SonNReal = %d) !!\n",
                    SonLv, SonPID0, amr->NPatchComma[SonLv][1] );

      if ( SonPID0%8 != 0 )
         Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 (%d) is not a multiple of 8 !!\n", SonLv, SonPID0 );
   }
#  endif


// 3. store and sort the 1D corner coordinates for all target son patches with LocalID == 0
   for (int t=0; t<NTargetSon0; t++)
   {
      SonPID0      = TargetSonPID0[t];
      Cr1D_Son0[t] = amr->patch[0][SonLv][SonPID0]->PaddedCr1D;
   }

   Mis_Heapsort( NTargetSon0, Cr1D_Son0, Cr1D_Son0_IdxTable );


// 4. matching
   Mis_Matching_int( amr->num[FaLv], amr->LB->PaddedCr1DList[FaLv], NTargetSon0, Cr1D_Son0, Match_Son0 );


// 5. construct father <-> son relation
   for (int t=0; t<NTargetSon0; t++)
   {
      SonPID0 = TargetSonPID0[ Cr1D_Son0_IdxTable[t] ];
      FaID    = Match_Son0[t];

      if ( FaID != -1 ) // father is found
      {
         FaPID = amr->LB->PaddedCr1DList_IdxTable[FaLv][FaID];

//       son -> father
         for (SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)   amr->patch[0][SonLv][SonPID]->father = FaPID;

//       father -> son
         amr->patch[0][FaLv][FaPID]->son = SonPID0;
      }

#     ifdef GAMER_DEBUG
      else // find no father (should NOT happen for any real patches)
         Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d found no father !!\n", SonLv, SonPID0 );
#     endif
   }


// 6. check results in debug mode
#  ifdef GAMER_DEBUG
   const int FaNNoFaBuf = amr->NPatchComma[FaLv][2];  // exclude father-buffer patches

// check 1 : all real patches must have fathers
   for (int t=0; t<NTargetSon0; t++)
   {
      SonPID0 = TargetSonPID0[t];

      for (SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
         if ( amr->patch[0][SonLv][SonPID]->father == -1 )
            Aux_Error( ERROR_INFO, "Check 1, SonLv %d: SonPID (%d) has no father !!\n", SonLv, SonPID );
   }

// check 2 : SonPIDs should always be real patches at SonLv
   for (FaPID=0; FaPID<FaNNoFaBuf; FaPID++)
   {
      SonPID = amr->patch[0][FaLv][FaPID]->son;

      if ( SonPID >= amr->NPatchComma[SonLv][1] )
         Aux_Error( ERROR_INFO, "Check 2, FaLv %d: FaPID (%d) -> son (%d) >= SonNReal (%d) !!\n",
                    FaLv, FaPID, SonPID, amr->NPatchComma[SonLv][1] );
   }

// check 3 : father's son = itself
   for (int t=0; t<NTargetSon0; t++)
   {
      SonPID0 = TargetSonPID0[t];

      for (SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
      {
         FaPID = amr->patch[0][SonLv][SonPID]->father;

         if ( amr->patch[0][FaLv][FaPID]->son != SonPID0 )
            Aux_Error( ERROR_INFO, "Check 3, SonLv %d: FaPID (%d)->son (%d) != SonPID (%d)->SonPID0 (%d) !!\n",
                       SonLv, FaPID, amr->patch[0][FaLv][FaPID]->son, SonPID, SonPID0 );
      }
   }

// check 4 : son's father = itself
   for (FaPID=0; FaPID<FaNNoFaBuf; FaPID++)
   {
      SonPID0 = amr->patch[0][FaLv][FaPID]->son;

      if ( SonPID0 >= 0 )
      {
         for (SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
         {
            if ( amr->patch[0][SonLv][SonPID]->father != FaPID )
               Aux_Error( ERROR_INFO, "Check 4, SonLv %d: SonPID (%d) -> father (%d) != FaPID (%d) !!\n",
                          SonLv, SonPID, amr->patch[0][SonLv][SonPID]->father, FaPID );
         }
      }
   }

// check 5 : double check the father <-> son relation
   int *Corner_Son, *Corner_Fa;

   for (int t=0; t<NTargetSon0; t++)
   {
      SonPID0    = TargetSonPID0[t];
      Corner_Son = amr->patch[0][SonLv][SonPID0]->corner;

      for (FaPID=0; FaPID<FaNNoFaBuf; FaPID++)
      {
         Corner_Fa = amr->patch[0][FaLv][FaPID]->corner;

         if (  Corner_Son[0] == Corner_Fa[0]  &&
               Corner_Son[1] == Corner_Fa[1]  &&
               Corner_Son[2] == Corner_Fa[2]     )
         {
            if ( amr->patch[0][SonLv][SonPID0]->father != FaPID )
               Aux_Error( ERROR_INFO, "Check 5.1, SonLv %d: SonPID0 (%d) -> father (%d) != FaPID (%d) !!\n",
                          SonLv, SonPID0, amr->patch[0][SonLv][SonPID0]->father, FaPID );

            if ( amr->patch[0][FaLv][FaPID]->son != SonPID0 )
               Aux_Error( ERROR_INFO, "Check 5.2, FaLv %d: FaPID (%d) -> son (%d) != SonPID0 (%d) !!\n",
                          FaLv, FaPID, amr->patch[0][FaLv][FaPID]->son, SonPID0 );
         }
      }
   }

// check 6 : father LB_Idx*8 == son LB_Idx - son LB_Idx%8
#  if ( LOAD_BALANCE == HILBERT )
   for (int t=0; t<NTargetSon0; t++)
   {
      SonPID0 = TargetSonPID0[t];

      for (SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
      {
         FaPID = amr->patch[0][SonLv][SonPID]->father;

         if ( amr->patch[0][FaLv][FaPID]->LB_Idx*8 != amr->patch[0][SonLv][SonPID]->LB_Idx - amr->patch[0][SonLv][SonPID]->LB_Idx%8 )
            Aux_Error( ERROR_INFO, "Check 6, SonLv %d: FaPID (%d)->LB_Idx (%ld) and SonPID (%d)->LB_Idx (%ld) are inconsistent !!\n",
                       SonLv, FaPID, amr->patch[0][FaLv][FaPID]->LB_Idx, SonPID, amr->patch[0][SonLv][SonPID]->LB_Idx );
      }
   }
#  endif
#  endif // #ifdef GAMER_DEBUG


// free memory
   delete [] Cr1D_Son0;
   delete [] Cr1D_Son0_IdxTable;
   delete [] Match_Son0;
   if ( SearchAllSon )  delete [] TargetSonPID0;

} // FUNCTION : LB_FindFather



#endif // #ifdef LOAD_BALANCE
