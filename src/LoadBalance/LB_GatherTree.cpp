#include "GAMER.h"

/*
Instructions for adding new patch_t members to GatherTree:
-> modify "patch.h"
-  add member to LB_GlobalPatch 
-  add lists to LB_LocalPatchExchangeList and LB_GlobalPatchExchangeList
-> modify "LB_GatherTree.cpp"
-  read new member in LB_FillLocalExchangeList
-  transfer member in LB_FillGlobalExchangeList 
-  write member to LB_GlobalPatch in LB_ConstructGlobalTree
*/


LB_PatchCount::LB_PatchCount() : NPatchAllLv(0), isInitialised(false) {
   NPatchAllRank = new int [MPI_NRank][NLEVEL];
   for (int lv=0; lv<NLEVEL; lv++)
   {
      GID_Offset[lv] = 0; 
      GID_LvStart[lv] = 0;
   } 
}

LB_PatchCount::~LB_PatchCount() {
   delete [] NPatchAllRank;
}

LB_LocalPatchExchangeList::LB_LocalPatchExchangeList() : isInitialised(false), LBIdxisInitialised(false) {
   // local lists for storing local tree structure
   for (int lv=0; lv<NLEVEL; lv++)
   {
      LBIdxList_Local        [lv] = new long [ amr->NPatchComma[lv][1] ];
      CrList_Local           [lv] = new int  [ amr->NPatchComma[lv][1] ][3];
      FaList_Local           [lv] = new int  [ amr->NPatchComma[lv][1] ];
      SonList_Local          [lv] = new int  [ amr->NPatchComma[lv][1] ];
      SibList_Local          [lv] = new int  [ amr->NPatchComma[lv][1] ][26];
#        ifdef PARTICLE
      NParList_Local         [lv] = new int  [ amr->NPatchComma[lv][1] ];
#        endif

      LBIdxList_Sort         [lv] = new long [ NPatchTotal[lv] ];
      LBIdxList_Sort_IdxTable[lv] = new int  [ NPatchTotal[lv] ];
   }

}

LB_LocalPatchExchangeList::~LB_LocalPatchExchangeList() { 
   for (int lv=0; lv<NLEVEL; lv++)
   {
      delete [] LBIdxList_Local[lv];
      delete []    CrList_Local[lv];
      delete []    FaList_Local[lv];
      delete []   SonList_Local[lv];
      delete []   SibList_Local[lv];
#        ifdef PARTICLE
      delete []  NParList_Local[lv];
#        endif

      delete [] LBIdxList_Sort [lv];
      delete [] LBIdxList_Sort_IdxTable[lv];
   }
}


// allocate memory and store pointers to lists with global patch information
LB_GlobalPatchExchangeList::LB_GlobalPatchExchangeList(LB_PatchCount& pc, int root) : isAllocated(false), isInitialised(false) {
// allocate lists for all ranks or for root rank
   if ( root < 0 || root == MPI_Rank) {
      LBIdxList_AllLv = new long [ pc.NPatchAllLv ];
      CrList_AllLv    = new int  [ pc.NPatchAllLv ][3];
      FaList_AllLv    = new int  [ pc.NPatchAllLv ];
      SonList_AllLv   = new int  [ pc.NPatchAllLv ];
      SibList_AllLv   = new int  [ pc.NPatchAllLv ][26];
#     ifdef PARTICLE
      NParList_AllLv  = new int  [ pc.NPatchAllLv ];
#     endif // # ifdef PARTICLE
      isAllocated     = true;
   } else {
      LBIdxList_AllLv = NULL; 
      CrList_AllLv    = NULL; 
      FaList_AllLv    = NULL; 
      SonList_AllLv   = NULL; 
      SibList_AllLv   = NULL; 
#     ifdef PARTICLE
      NParList_AllLv  = NULL; 
#     endif // # ifdef PARTICLE
   }
}

LB_GlobalPatchExchangeList::~LB_GlobalPatchExchangeList() {
   if ( isAllocated ) {
      delete [] LBIdxList_AllLv;
      delete []    CrList_AllLv;
      delete []    FaList_AllLv;
      delete []   SonList_AllLv;
      delete []   SibList_AllLv;
#     ifdef PARTICLE
      delete []  NParList_AllLv;
#     endif // #ifdef PARTICLE
   }
}

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GetPID
// Description :  Convert GID to local PID
//
// Note        :  1. Calculate PID and level from PID.
//
// Parameter   :  GID      : GID to convert
//                level    : fill in level of GID.
//                PID      : fill in PID of GID
//                rank     : fill in rank of GID
//
// Return      :  *level, *PID
//-------------------------------------------------------------------------------------------------------
void LB_GetPID(long GID, int& level, int& PID, LB_PatchCount& pc) {
#  ifdef GAMER_DEBUG
   if ( GID < 0  ||  GID >= pc.NPatchAllLv )  
      Aux_Error( ERROR_INFO, "incorrect GID %ld (max = %ld) !!\n", GID, pc.NPatchAllLv-1 );
#  endif

   level = 0;

   for(int lv = 1; lv < NLEVEL; lv++) {
      if ( GID < pc.GID_Offset[lv] )      
        break;
      level = lv;
   }

   PID = GID - pc.GID_Offset[level];
}



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_AllgatherPatchCount
// Description :  Gather the number of patches at different MPI ranks and set the corresponding GID offset
//
// Note        :  - store data in LB_PatchCount& pc
//
// Parameter   :  pc   : reference to LB_PatchCount object
//-------------------------------------------------------------------------------------------------------
void LB_AllgatherPatchCount(LB_PatchCount& pc) {

   pc.NPatchLocalLv = 0; 
   pc.NPatchAllLv   = 0;

   for (int lv=0; lv<NLEVEL; lv++)  {
      pc.NPatchLocal[lv] = amr->NPatchComma[lv][1];
      pc.NPatchLocalLv += pc.NPatchLocal[lv];
   }

   MPI_Allgather( pc.NPatchLocal, NLEVEL, MPI_INT, pc.NPatchAllRank[0], NLEVEL, MPI_INT, MPI_COMM_WORLD );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      pc.GID_Offset[lv] = 0;

      for (int r=0; r<MPI_Rank; r++)      pc.GID_Offset[lv] += pc.NPatchAllRank[r][lv];

      for (int FaLv=0; FaLv<lv; FaLv++)   pc.GID_Offset[lv] += NPatchTotal[FaLv];

      pc.NPatchAllLv += NPatchTotal[lv];

      pc.GID_LvStart[lv] = ( lv == 0 ) ? 0 : pc.GID_LvStart[lv-1] + NPatchTotal[lv-1];
   }

   pc.isInitialised = true; 
}


//-------------------------------------------------------------------------------------------------------
// Function    :  LB_AllgatherLBIdx
// Description :  Collect and sort LBIdx from all ranks
//
// Note        :  - pc requires initialisation by calling LB_AllgatherPatchCount
//
// Parameter   :  pc   : reference to LB_PatchCount object
//             :  lel  : reference to LB_LocalPatchExchangeList
//             :  gel  : reference to LB_GlobalPatchExchangeList
//             :  root : root MPI rank, -1 for all ranks
//-------------------------------------------------------------------------------------------------------
void LB_AllgatherLBIdx(LB_PatchCount& pc, LB_LocalPatchExchangeList& lel, LB_GlobalPatchExchangeList* gel)  {

#  ifdef GAMER_DEBUG
   if ( !pc.isInitialised ) 
      Aux_Error( ERROR_INFO, "call LB_AllgatherLBIdx without initialising LB_PatchCount object !!\n");
#  endif 

   int RecvCount_LBIdx[MPI_NRank], RecvDisp_LBIdx[MPI_NRank];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvCount_LBIdx[r] = pc.NPatchAllRank[r][lv];
         RecvDisp_LBIdx [r] = ( r == 0 ) ? 0 : RecvDisp_LBIdx[r-1] + RecvCount_LBIdx[r-1];
      }

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         lel.LBIdxList_Local[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

//    all ranks need to get LBIdxList_Sort since we will use it to calculate GID
      MPI_Allgatherv( lel.LBIdxList_Local[lv], amr->NPatchComma[lv][1], MPI_LONG,
                      lel.LBIdxList_Sort[lv], RecvCount_LBIdx, RecvDisp_LBIdx, MPI_LONG,
                      MPI_COMM_WORLD );
   } // for (int lv=0; lv<NLEVEL; lv++)

// store in the AllLv array BEFORE sorting
   if ( gel != NULL ) {
      if ( gel->isAllocated ) {
         long MyGID = 0;

         for (int lv=0; lv<NLEVEL; lv++)
         for (int PID=0; PID<NPatchTotal[lv]; PID++)
            gel->LBIdxList_AllLv[ MyGID++ ] = lel.LBIdxList_Sort[lv][PID];
      }
   }

// sort list and get the corresponding index table (for calculating GID later)
   for (int lv=0; lv<NLEVEL; lv++)
      Mis_Heapsort( NPatchTotal[lv], lel.LBIdxList_Sort[lv], lel.LBIdxList_Sort_IdxTable[lv] );

   lel.LBIdxisInitialised = true; 
}

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_FillLocalPatchExchangeList
// Description :  Fill local exchange list by reading amr->patch structure on local MPI rank
//
// Note        :  - pc requires initialisation by calling LB_AllgatherPatchCount
//                - lel requires initialisation by calling LB_AllgatherLBIdx
//
// Parameter   :  pc   : reference to LB_PatchCount object
//             :  lel  : reference to LB_LocalPatchExchangeList
//-------------------------------------------------------------------------------------------------------
void LB_FillLocalPatchExchangeList(LB_PatchCount& pc, LB_LocalPatchExchangeList& lel) {

#  ifdef GAMER_DEBUG
   if ( !pc.isInitialised ) 
      Aux_Error( ERROR_INFO, "call LB_FillLocalExchangeList without initialising LB_PatchCount object !!\n");
   if ( !lel.LBIdxisInitialised ) 
      Aux_Error( ERROR_INFO, "call LB_FillLocalExchangeList without initialising load balancing id lists object !!\n");
#  endif 

// temporary variables
   int   MyGID, FaPID, FaGID, FaLv, SonPID, SonGID, SonLv, SibPID, SibGID, MatchIdx;
   long  FaLBIdx, SonLBIdx, SibLBIdx;
   int  *SonCr=NULL, *SibCr=NULL;

// 4. store the local tree
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       4-1. LBIdx (set already)
//       lel.LBIdxList_Local[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;


//       4-2. corner
         for (int d=0; d<3; d++)
         lel.CrList_Local[lv][PID][d] = amr->patch[0][lv][PID]->corner[d];


//       4-3. father GID
         FaPID = amr->patch[0][lv][PID]->father;
         FaLv  = lv - 1;

//       no father (only possible for the root patches)
         if ( FaPID < 0 )
         {
#           ifdef GAMER_DEBUG
            if ( lv != 0 )       Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d < 0 !!\n", lv, PID, FaPID );
            if ( FaPID != -1 )   Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d < 0 but != -1 !!\n", lv, PID, FaPID );
#           endif

            FaGID = FaPID;
         }

//       father patch is a real patch
         else if ( FaPID < amr->NPatchComma[FaLv][1] )
            FaGID = FaPID + pc.GID_Offset[FaLv];

//       father patch is a buffer patch (only possible in LOAD_BALANCE)
         else // (FaPID >= amr->NPatchComma[FaLv][1] )
         {
#           ifdef GAMER_DEBUG
#           ifndef LOAD_BALANCE
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d >= NRealFaPatch %d (only possible in LOAD_BALANCE) !!\n",
                       lv, PID, FaPID, amr->NPatchComma[FaLv][1] );
#           endif

            if ( FaPID >= amr->num[FaLv] )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d >= total number of patches %d !!\n",
                       lv, PID, FaPID, amr->num[FaLv] );
#           endif // GAMER_DEBUG

            FaLBIdx = amr->patch[0][FaLv][FaPID]->LB_Idx;

            Mis_Matching_int( NPatchTotal[FaLv], lel.LBIdxList_Sort[FaLv], 1, &FaLBIdx, &MatchIdx );

#           ifdef GAMER_DEBUG
            if ( MatchIdx < 0 )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d, FaLBIdx %ld, couldn't find a matching patch !!\n",
                       lv, PID, FaPID, FaLBIdx );
#           endif

            FaGID = lel.LBIdxList_Sort_IdxTable[FaLv][MatchIdx] + pc.GID_LvStart[FaLv];
         } // if ( FaPID >= amr->NPatchComma[FaLv][1] )

         lel.FaList_Local[lv][PID] = FaGID;


//       4-4. son GID
         SonPID = amr->patch[0][lv][PID]->son;
         SonLv  = lv + 1;

//       no son (must check this first since SonLv may be out of range --> == NLEVEL)
         if      ( SonPID == -1 )
            SonGID = SonPID;

//       son patch is a real patch at home
         else if ( SonPID >= 0  &&  SonPID < amr->NPatchComma[SonLv][1] )
            SonGID = SonPID + pc.GID_Offset[SonLv];

//       son patch lives abroad (only possible in LOAD_BALANCE)
         else if ( SonPID < -1 )
         {
#           ifdef GAMER_DEBUG
#           ifdef LOAD_BALANCE
            const int SonRank = SON_OFFSET_LB - SonPID;
            if ( SonRank < 0  ||  SonRank == MPI_Rank  ||  SonRank >= MPI_NRank )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d, incorrect SonRank %d (MyRank %d, NRank %d) !!\n",
                       lv, PID, SonPID, SonRank, MPI_Rank, MPI_NRank );
#           else
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d < -1 (only possible in LOAD_BALANCE) !!\n",
                       lv, PID, SonPID );
#           endif // LOAD_BALANCE
#           endif // GAMER_DEBUG

//          get the SonGID by "father corner = son corner -> son LB_Idx -> son GID"
//          --> for Hilbert curve we have "SonLBIdx-SonLBIdx%8 = 8*MyLBIdx"
            SonCr    = amr->patch[0][lv][PID]->corner;
            SonLBIdx = LB_Corner2Index( SonLv, SonCr, CHECK_ON );

#           if ( defined GAMER_DEBUG  &&  LOAD_BALANCE == HILBERT )
            if ( SonLBIdx - SonLBIdx%8 != 8*amr->patch[0][lv][PID]->LB_Idx )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d, SonCr (%d,%d,%d), incorret SonLBIdx %ld, (MyLBIdx %ld) !!\n",
                       lv, PID, SonPID, SonCr[0], SonCr[1], SonCr[2], SonLBIdx, amr->patch[0][lv][PID]->LB_Idx );
#           endif

            Mis_Matching_int( NPatchTotal[SonLv], lel.LBIdxList_Sort[SonLv], 1, &SonLBIdx, &MatchIdx );

#           ifdef GAMER_DEBUG
            if ( MatchIdx < 0 )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d, SonLBIdx %ld, couldn't find a matching patch !!\n",
                       lv, PID, SonPID, SonLBIdx );
#           endif

            SonGID = lel.LBIdxList_Sort_IdxTable[SonLv][MatchIdx] + pc.GID_LvStart[SonLv];
         } // else if ( SonPID < -1 )

//       son patch is a buffer patch (SonPID >= amr->NPatchComma[SonLv][1]) --> impossible
         else // ( SonPID >= amr->NPatchComma[SonLv][1] )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d is a buffer patch (NRealSonPatch %d) !!\n",
                       lv, PID, SonPID, amr->NPatchComma[SonLv][1] );

         lel.SonList_Local[lv][PID] = SonGID;


//       4-5. sibling GID
         for (int s=0; s<26; s++)
         {
            SibPID = amr->patch[0][lv][PID]->sibling[s];

//          no sibling (SibPID can be either -1 or SIB_OFFSET_NONPERIODIC-BoundaryDirection)
            if      ( SibPID < 0 )
               SibGID = SibPID;

//          sibling patch is a real patch
            else if ( SibPID < amr->NPatchComma[lv][1] )
               SibGID = SibPID + pc.GID_Offset[lv];

//          sibling patch is a buffer patch (which may lie outside the simulation domain)
            else
            {
#              ifdef GAMER_DEBUG
               if ( SibPID >= amr->num[lv] )
               Aux_Error( ERROR_INFO, "Lv %d, PID %d, SibPID %d >= total number of patches %d !!\n",
                          lv, PID, SibPID, amr->num[lv] );
#              endif

//             get the SibGID by "sibling corner -> sibling LB_Idx -> sibling GID"
               SibCr    = amr->patch[0][lv][SibPID]->corner;
               SibLBIdx = LB_Corner2Index( lv, SibCr, CHECK_OFF );   // periodicity has been assumed here

               Mis_Matching_int( NPatchTotal[lv], lel.LBIdxList_Sort[lv], 1, &SibLBIdx, &MatchIdx );

#              ifdef GAMER_DEBUG
               if ( MatchIdx < 0 )
               Aux_Error( ERROR_INFO, "Lv %d, PID %d, SibPID %d, SibLBIdx %ld, couldn't find a matching patch !!\n",
                          lv, PID, SibPID, SibLBIdx );
#              endif

               SibGID = lel.LBIdxList_Sort_IdxTable[lv][MatchIdx] + pc.GID_LvStart[lv];
            } // if ( SibPID >= amr->NPatchComma[lv][1] )

            lel.SibList_Local[lv][PID][s] = SibGID;

         } // for (int s=0; s<26; s++)


#        ifdef PARTICLE
//       4-6. NPar
         lel.NParList_Local[lv][PID] = amr->patch[0][lv][PID]->NPar;
#        endif

      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)

   lel.isInitialised = true; 
}

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_FillGlobalPatchExchangeList
// Description :  Fill global patch exchange lists by exchanging local patch list data between ranks
//
// Note        :  - pc and lel need to be initialised by calling LB_AllgatherPatchCount and LB_FillLocalPatchExchangeList beforehand
//                - global patch data is written to LB_GlobalPatchExchangeList& gel if root = gel.root or root = -1
//                - pass root = -1 to exchange local patch list data from all ranks to all ranks
//
// Parameter   :  pc   : reference to LB_PatchCount object
//             :  lel  : reference to LB_LocalPatchExchangeList
//             :  gel  : reference to LB_GlobalPatchExchangeList
//             :  root : root MPI rank, -1 for all ranks
//-------------------------------------------------------------------------------------------------------
void LB_FillGlobalPatchExchangeList(LB_PatchCount& pc, LB_LocalPatchExchangeList& lel, LB_GlobalPatchExchangeList& gel, int root) {
#  ifdef GAMER_DEBUG
   if ( !pc.isInitialised ) 
      Aux_Error( ERROR_INFO, "call LB_FillGlobalPatchExchangeList without initialising LB_PatchCount object !!\n");
   if ( !lel.isInitialised ) 
      Aux_Error( ERROR_INFO, "call LB_FillGlobalPatchExchangeList without initialising LB_LocalPatchExchangeList object !!\n");
#  endif 

// sending and receiving lists for MPI communication 
   int   RecvCount_Cr[MPI_NRank],    RecvDisp_Cr[MPI_NRank];
   int   RecvCount_Fa[MPI_NRank],    RecvDisp_Fa[MPI_NRank];
   int   RecvCount_Son[MPI_NRank],   RecvDisp_Son[MPI_NRank];
   int   RecvCount_Sib[MPI_NRank],   RecvDisp_Sib[MPI_NRank];
#  ifdef PARTICLE
   int   RecvCount_NPar[MPI_NRank],  RecvDisp_NPar[MPI_NRank];
#  endif // #  ifdef PARTICLE


// 5. gather data from all ranks
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvCount_Fa  [r] = pc.NPatchAllRank[r][lv];
         RecvCount_Son [r] = RecvCount_Fa[r];
         RecvCount_Sib [r] = RecvCount_Fa[r]*26;
         RecvCount_Cr  [r] = RecvCount_Fa[r]*3;
#        ifdef PARTICLE
         RecvCount_NPar[r] = RecvCount_Fa[r];
#        endif

         RecvDisp_Fa   [r] = ( r == 0 ) ? 0 : RecvDisp_Fa[r-1] + RecvCount_Fa[r-1];
         RecvDisp_Son  [r] = RecvDisp_Fa[r];
         RecvDisp_Sib  [r] = RecvDisp_Fa[r]*26;
         RecvDisp_Cr   [r] = RecvDisp_Fa[r]*3;
#        ifdef PARTICLE
         RecvDisp_NPar [r] = RecvDisp_Fa[r];
#        endif

      }

//    note that we collect data at one level at a time
      if (root < 0) {
         MPI_Allgatherv( lel.FaList_Local[lv],     amr->NPatchComma[lv][1],    MPI_INT,
                      gel.FaList_AllLv+pc.GID_LvStart[lv],       RecvCount_Fa,   RecvDisp_Fa,   MPI_INT, MPI_COMM_WORLD );

         MPI_Allgatherv( lel.SonList_Local[lv],    amr->NPatchComma[lv][1],    MPI_INT,
                      gel.SonList_AllLv+pc.GID_LvStart[lv],      RecvCount_Son,  RecvDisp_Son,  MPI_INT, MPI_COMM_WORLD );

         MPI_Allgatherv( lel.SibList_Local[lv][0], amr->NPatchComma[lv][1]*26, MPI_INT,
                      (gel.SibList_AllLv+pc.GID_LvStart[lv])[0], RecvCount_Sib,  RecvDisp_Sib,  MPI_INT, MPI_COMM_WORLD );

         MPI_Allgatherv( lel.CrList_Local[lv][0],  amr->NPatchComma[lv][1]*3,  MPI_INT,
                      (gel.CrList_AllLv+pc.GID_LvStart[lv])[0],  RecvCount_Cr,   RecvDisp_Cr,   MPI_INT, MPI_COMM_WORLD );

#        ifdef PARTICLE
         MPI_Allgatherv( lel.NParList_Local[lv],    amr->NPatchComma[lv][1],    MPI_INT,
                      gel.NParList_AllLv+pc.GID_LvStart[lv],     RecvCount_NPar, RecvDisp_NPar, MPI_INT, MPI_COMM_WORLD );
#        endif
      } else {
         MPI_Gatherv( lel.FaList_Local[lv],     amr->NPatchComma[lv][1],    MPI_INT,
                      gel.FaList_AllLv+pc.GID_LvStart[lv],       RecvCount_Fa,   RecvDisp_Fa,   MPI_INT, root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.SonList_Local[lv],    amr->NPatchComma[lv][1],    MPI_INT,
                      gel.SonList_AllLv+pc.GID_LvStart[lv],      RecvCount_Son,  RecvDisp_Son,  MPI_INT, root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.SibList_Local[lv][0], amr->NPatchComma[lv][1]*26, MPI_INT,
                      (gel.SibList_AllLv+pc.GID_LvStart[lv])[0], RecvCount_Sib,  RecvDisp_Sib,  MPI_INT, root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.CrList_Local[lv][0],  amr->NPatchComma[lv][1]*3,  MPI_INT,
                      (gel.CrList_AllLv+pc.GID_LvStart[lv])[0],  RecvCount_Cr,   RecvDisp_Cr,   MPI_INT, root, MPI_COMM_WORLD );

#        ifdef PARTICLE
         MPI_Gatherv( lel.NParList_Local[lv],    amr->NPatchComma[lv][1],    MPI_INT,
                      gel.NParList_AllLv+pc.GID_LvStart[lv],     RecvCount_NPar, RecvDisp_NPar, MPI_INT, root, MPI_COMM_WORLD );
#        endif
      }

   } // for (int lv=0; lv<NLEVEL; lv++)

   gel.isInitialised = true; 
}

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_ConstructGlobalTree
// Description :  Gather global tree structure as vector indexed by GIDs to root rank
//
// Note        :  - Store global tree AMR structure gathered from all ranks in vector 
//                - Example usage: Print level of MPI rank 0's patch's son
//                      LB_PatchCount pc;
//                      std::vector<LB_GlobalPatch> t = LB_GatherTree(pc, 0);
//                      if (t[GID].son != -1) printf(t[t[GID].son].level);
//
// Parameter   :  pc   : reference to LB_PatchCount object
//             :  gel  : reference to LB_GlobalPatchExchangeList that needs to be initialised by calling LB_FillGlobalPatchExchangeList
//             :  root : root MPI rank that receives global list, -1 for all ranks
//-------------------------------------------------------------------------------------------------------
LB_GlobalPatch* LB_ConstructGlobalTree(LB_PatchCount& pc, LB_GlobalPatchExchangeList& gel, int root) {
   LB_GlobalPatch* global_tree;
   if ( root >= 0 && root != MPI_Rank )
      return global_tree; 

#  ifdef GAMER_DEBUG
   if ( !pc.isInitialised ) 
      Aux_Error( ERROR_INFO, "call LB_ConstructGlobalTree without initialising LB_PatchCount object !!\n");
   if ( !gel.isInitialised ) 
      Aux_Error( ERROR_INFO, "call LB_ConstructGlobalTree without initialising LB_GlobalPatchExchangeListobject !!\n");
#  endif 

   global_tree = new LB_GlobalPatch[pc.NPatchAllLv];

   long MyGID = 0;
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int i = 0; i < NPatchTotal[lv]; ++i) {
         global_tree[MyGID].level      = lv; 
         global_tree[MyGID].father     = gel.FaList_AllLv   [MyGID];
         global_tree[MyGID].son        = gel.SonList_AllLv  [MyGID];
         for (int s=0; s<26; s++) 
         global_tree[MyGID].sibling[s] = gel.SibList_AllLv  [MyGID][s];
         for (int c=0; c<3 ; c++) 
         global_tree[MyGID].corner[c]  = gel.CrList_AllLv   [MyGID][c];
#        ifdef PARTICLE
         global_tree[MyGID].NPar       = gel.NParList_AllLv [MyGID];
#        endif //# ifdef PARTICLE
         global_tree[MyGID].LB_Idx     = gel.LBIdxList_AllLv[MyGID];
         MyGID += 1; 
      }
   }

   return global_tree;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GatherTree
// Description :  Gather global tree structure as vector indexed by GIDs to root rank
//
// Note        :  - Store global tree AMR structure gathered from all ranks in vector 
//                - Example usage: Print level of MPI rank 0's patch's son
//                      LB_PatchCount pc;
//                      std::vector<LB_GlobalPatch> t = LB_GatherTree(pc, 0);
//                      if (t[GID].son != -1) printf(t[t[GID].son].level);
//
// Parameter   :  pc  : reference to LB_PatchCount object that is filled with the GID offsets for 
//                      the local rank as well as total number of patches
//-------------------------------------------------------------------------------------------------------

LB_GlobalPatch* LB_GatherTree(LB_PatchCount& pc, int root) {
// get patch counts per level and per rank from all ranks
   LB_AllgatherPatchCount(pc); 

// set up local and global exchange lists
   LB_LocalPatchExchangeList  lel;
   LB_GlobalPatchExchangeList gel(pc, root);

// exchange load balance id's between all ranks
   LB_AllgatherLBIdx(pc, lel, &gel);

// fill local patch lists with information from patches
   LB_FillLocalPatchExchangeList(pc, lel);

// exchange local patch information with other ranks
   LB_FillGlobalPatchExchangeList(pc, lel, gel, root);

// construct and return vector with global tree information
   return LB_ConstructGlobalTree(pc, gel, root);
}

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_AllgatherTree
// Description :  Gather global tree structure as vector indexed by GIDs to all ranks
//
// Note        :  - see LB_GatherTrees
//
// Parameter   :  pc  : reference to LB_PatchCount object that is filled with the GID offsets for the local rank as well as total number of patches
//-------------------------------------------------------------------------------------------------------
LB_GlobalPatch* LB_AllgatherTree(LB_PatchCount& pc) {
// call LB_GatherTree with root = -1 to synchronise data between all ranks
   return LB_GatherTree(pc, -1);
}