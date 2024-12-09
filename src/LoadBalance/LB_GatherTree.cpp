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



LB_PatchCount::LB_PatchCount() : NPatchAllLv(0), NPatchLocalAllLv(0), isInitialised(false) {

   NPatchAllRank = new int [MPI_NRank][NLEVEL];
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_Rank; r++)   NPatchAllRank[r][lv] = 0;
      NPatchLocal[lv] = 0;
      GID_Offset [lv] = 0;
      GID_LvStart[lv] = 0;
   }

} // FUNCTION : LB_PatchCount



LB_PatchCount::~LB_PatchCount() {

   delete [] NPatchAllRank;

} // FUNCTION : ~LB_PatchCount



LB_LocalPatchExchangeList::LB_LocalPatchExchangeList() : isInitialised(false), LBIdxisInitialised(false) {

// local lists for storing local tree structure
   for (int lv=0; lv<NLEVEL; lv++)
   {
      LBIdxList_Local        [lv] = new long   [ amr->NPatchComma[lv][1] ];
      CrList_Local           [lv] = new int    [ amr->NPatchComma[lv][1] ][3];
      FaList_Local           [lv] = new int    [ amr->NPatchComma[lv][1] ];
      SonList_Local          [lv] = new int    [ amr->NPatchComma[lv][1] ];
      SibList_Local          [lv] = new int    [ amr->NPatchComma[lv][1] ][26];
      EdgeLList_Local        [lv] = new double [ amr->NPatchComma[lv][1] ][3];
      EdgeRList_Local        [lv] = new double [ amr->NPatchComma[lv][1] ][3];
      PaddedCr1DList_Local   [lv] = new ulong  [ amr->NPatchComma[lv][1] ];
      MPI_RankList_Local     [lv] = new int    [ amr->NPatchComma[lv][1] ];
#     ifdef PARTICLE
      NParList_Local         [lv] = new int    [ amr->NPatchComma[lv][1] ];
#     endif

      LBIdxList_Sort         [lv] = new long   [ NPatchTotal[lv] ];
      LBIdxList_Sort_IdxTable[lv] = new int    [ NPatchTotal[lv] ];
   }

} // FUNCTION : LB_LocalPatchExchangeList



LB_LocalPatchExchangeList::~LB_LocalPatchExchangeList() {

   for (int lv=0; lv<NLEVEL; lv++)
   {
      delete []         LBIdxList_Local[lv];
      delete []            CrList_Local[lv];
      delete []            FaList_Local[lv];
      delete []           SonList_Local[lv];
      delete []           SibList_Local[lv];
      delete []         EdgeLList_Local[lv];
      delete []         EdgeRList_Local[lv];
      delete []    PaddedCr1DList_Local[lv];
      delete []      MPI_RankList_Local[lv];
#     ifdef PARTICLE
      delete []          NParList_Local[lv];
#     endif
      delete [] LBIdxList_Sort         [lv];
      delete [] LBIdxList_Sort_IdxTable[lv];
   }

} // FUNCTION : ~LB_LocalPatchExchangeList



// allocate memory and store pointers to lists with global patch information
LB_GlobalPatchExchangeList::LB_GlobalPatchExchangeList( LB_PatchCount& pc, int root ) : isAllocated(false), isInitialised(false) {

#  ifdef GAMER_DEBUG
   if ( !pc.isInitialised )
      Aux_Error( ERROR_INFO, "create object of type LB_GlobalPatchExchangeList without initialising LB_PatchCount object !!\n");
#  endif

// allocate lists for all ranks or for root rank
   if ( root < 0  ||  root == MPI_Rank) {
      LBIdxList_AllLv      = new long   [ pc.NPatchAllLv ];
      CrList_AllLv         = new int    [ pc.NPatchAllLv ][3];
      FaList_AllLv         = new int    [ pc.NPatchAllLv ];
      SonList_AllLv        = new int    [ pc.NPatchAllLv ];
      SibList_AllLv        = new int    [ pc.NPatchAllLv ][26];
      EdgeLList_AllLv      = new double [ pc.NPatchAllLv ][3];
      EdgeRList_AllLv      = new double [ pc.NPatchAllLv ][3];
      PaddedCr1DList_AllLv = new ulong  [ pc.NPatchAllLv ];
      MPI_RankList_AllLv   = new int    [ pc.NPatchAllLv ];
#     ifdef PARTICLE
      NParList_AllLv       = new int    [ pc.NPatchAllLv ];
#     endif

//    set allocation flag
      isAllocated          = true;
   } else {
      LBIdxList_AllLv      = NULL;
      CrList_AllLv         = NULL;
      FaList_AllLv         = NULL;
      SonList_AllLv        = NULL;
      SibList_AllLv        = NULL;
      EdgeLList_AllLv      = NULL;
      EdgeRList_AllLv      = NULL;
      PaddedCr1DList_AllLv = NULL;
      MPI_RankList_AllLv   = NULL;
#     ifdef PARTICLE
      NParList_AllLv       = NULL;
#     endif
   }

} // FUNCTION : LB_GlobalPatchExchangeList



LB_GlobalPatchExchangeList::~LB_GlobalPatchExchangeList() {

   if ( isAllocated ) {
      delete []      LBIdxList_AllLv;
      delete []         CrList_AllLv;
      delete []         FaList_AllLv;
      delete []        SonList_AllLv;
      delete []        SibList_AllLv;
      delete []      EdgeLList_AllLv;
      delete []      EdgeRList_AllLv;
      delete [] PaddedCr1DList_AllLv;
      delete []   MPI_RankList_AllLv;
#     ifdef PARTICLE
      delete []       NParList_AllLv;
#     endif
   }

} // FUNCTION : ~LB_GlobalPatchExchangeList



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GetPID
// Description :  Convert GID to local PID
//
// Note        :  - Calculate PID and level from GID.
//
// Parameter   :  GID        : GID to convert
//             :  level      : Reference to integer where level corresponding to GID is stored
//             :  PID        : Reference to integer where PID corrsponding to GID is stored
//             :  GID_Offset : Pointer to table with GID offsets on rank; array of length NLEVEL
//-------------------------------------------------------------------------------------------------------
void LB_GetPID( const int GID, int& level, int& PID, int* GID_Offset ) {

#   ifdef GAMER_DEBUG
    long NPatchAllLv = 0;
    for (int lv=0; lv<NLEVEL; lv++)    NPatchAllLv += NPatchTotal[lv];
    if ( GID < 0  ||  GID >= NPatchAllLv )   Aux_Error( ERROR_INFO, "incorrect gid %ld (max = %ld) !!\n", GID, NPatchAllLv-1 );
#   endif

   level = 0;

   for(int lv=1; lv<NLEVEL; lv++) {
      if ( GID < GID_Offset[lv] )
        break;
      level = lv;
   }

   PID = GID - GID_Offset[level];

} // FUNCTION : LB_GetPID



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_AllgatherPatchCount
// Description :  Gather the number of patches at different MPI ranks and set the corresponding GID offset
//
// Note        :  - Store data in LB_PatchCount& pc
//
// Parameter   :  pc : Reference to LB_PatchCount object
//-------------------------------------------------------------------------------------------------------
void LB_AllgatherPatchCount( LB_PatchCount& pc ) {

   pc.NPatchLocalAllLv = 0;
   pc.NPatchAllLv      = 0;

   for (int lv=0; lv<NLEVEL; lv++)  {
      pc.NPatchLocal[lv] = amr->NPatchComma[lv][1];
      pc.NPatchLocalAllLv += pc.NPatchLocal[lv];
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

} // FUNCTION : LB_AllgatherPatchCount



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_AllgatherLBIdx
// Description :  Collect and sort LBIdx from all ranks
//
// Note        :  - pc requires initialisation by calling LB_AllgatherPatchCount
//
// Parameter   :  pc   : Reference to LB_PatchCount object
//             :  lel  : Reference to LB_LocalPatchExchangeList
//             :  gel  : Reference to LB_GlobalPatchExchangeList
//             :  root : Root MPI rank, -1 for all ranks
//-------------------------------------------------------------------------------------------------------
void LB_AllgatherLBIdx( LB_PatchCount& pc, LB_LocalPatchExchangeList& lel, LB_GlobalPatchExchangeList* gel )  {

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
         int MyGID = 0;

         for (int lv=0; lv<NLEVEL; lv++)
         for (int PID=0; PID<NPatchTotal[lv]; PID++)
            gel->LBIdxList_AllLv[ MyGID++ ] = lel.LBIdxList_Sort[lv][PID];
      }
   }

// sort list and get the corresponding index table (for calculating GID later)
   for (int lv=0; lv<NLEVEL; lv++)
      Mis_Heapsort( NPatchTotal[lv], lel.LBIdxList_Sort[lv], lel.LBIdxList_Sort_IdxTable[lv] );

   lel.LBIdxisInitialised = true;

} // FUNCTION : LB_AllgatherLBIdx



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_FillLocalPatchExchangeList
// Description :  Fill local exchange list by reading amr->patch structure on local MPI rank
//
// Note        :  - pc requires initialisation by calling LB_AllgatherPatchCount
//                - lel requires initialisation by calling LB_AllgatherLBIdx
//
// Parameter   :  pc  : Reference to LB_PatchCount object
//             :  lel : Reference to LB_LocalPatchExchangeList
//-------------------------------------------------------------------------------------------------------
void LB_FillLocalPatchExchangeList( LB_PatchCount& pc, LB_LocalPatchExchangeList& lel ) {

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

// store the local tree
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       1. LBIdx (set already)
//       lel.LBIdxList_Local[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;


//       2. corner
         for (int d=0; d<3; d++)
         lel.CrList_Local[lv][PID][d] = amr->patch[0][lv][PID]->corner[d];


//       3. father GID
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


//       4. son GID
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


//       5. sibling GID
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

//       6. left edge
         for (int d=0; d<3; d++)
         lel.EdgeLList_Local[lv][PID][d] = amr->patch[0][lv][PID]->EdgeL[d];

//       7. right edge
         for (int d=0; d<3; d++)
         lel.EdgeRList_Local[lv][PID][d] = amr->patch[0][lv][PID]->EdgeR[d];

//       8. PaddedCr1D
         lel.PaddedCr1DList_Local[lv][PID] = amr->patch[0][lv][PID]->PaddedCr1D;

//       9. MPI Rank
         lel.MPI_RankList_Local[lv][PID] = MPI_Rank;

#        ifdef PARTICLE
//       10. NPar
         lel.NParList_Local[lv][PID] = amr->patch[0][lv][PID]->NPar;
#        endif

      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)

   lel.isInitialised = true;

} // FUNCTION : LB_FillLocalPatchExchangeList



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_FillGlobalPatchExchangeList
// Description :  Fill global patch exchange lists by exchanging local patch list data between ranks
//
// Note        :  - pc and lel need to be initialised by calling LB_AllgatherPatchCount and LB_FillLocalPatchExchangeList beforehand
//                - Global patch data is written to LB_GlobalPatchExchangeList& gel if root = gel.root or root = -1
//                - Pass root = -1 to exchange local patch list data from all ranks to all ranks
//
// Parameter   :  pc   : Reference to LB_PatchCount object
//             :  lel  : Reference to LB_LocalPatchExchangeList
//             :  gel  : Reference to LB_GlobalPatchExchangeList
//             :  root : Root MPI rank, -1 for all ranks
//-------------------------------------------------------------------------------------------------------
void LB_FillGlobalPatchExchangeList( LB_PatchCount& pc, LB_LocalPatchExchangeList& lel, LB_GlobalPatchExchangeList& gel, int root ) {

#  ifdef GAMER_DEBUG
   if ( !pc.isInitialised )
      Aux_Error( ERROR_INFO, "call LB_FillGlobalPatchExchangeList without initialising LB_PatchCount object !!\n");
   if ( !lel.isInitialised )
      Aux_Error( ERROR_INFO, "call LB_FillGlobalPatchExchangeList without initialising LB_LocalPatchExchangeList object !!\n");
#  endif

// sending and receiving lists for MPI communication
   int RecvCount_Cr        [MPI_NRank], RecvDisp_Cr        [MPI_NRank];
   int RecvCount_Fa        [MPI_NRank], RecvDisp_Fa        [MPI_NRank];
   int RecvCount_Son       [MPI_NRank], RecvDisp_Son       [MPI_NRank];
   int RecvCount_Sib       [MPI_NRank], RecvDisp_Sib       [MPI_NRank];
   int RecvCount_EdgeL     [MPI_NRank], RecvDisp_EdgeL     [MPI_NRank];
   int RecvCount_EdgeR     [MPI_NRank], RecvDisp_EdgeR     [MPI_NRank];
   int RecvCount_PaddedCr1D[MPI_NRank], RecvDisp_PaddedCr1D[MPI_NRank];
   int RecvCount_MPI_Rank  [MPI_NRank], RecvDisp_MPI_Rank  [MPI_NRank];
#  ifdef PARTICLE
   int RecvCount_NPar      [MPI_NRank], RecvDisp_NPar      [MPI_NRank];
#  endif


// gather data from all ranks
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvCount_Fa        [r] = pc.NPatchAllRank[r][lv];
         RecvCount_Son       [r] = RecvCount_Fa[r];
         RecvCount_Sib       [r] = RecvCount_Fa[r]*26;
         RecvCount_Cr        [r] = RecvCount_Fa[r]*3;
         RecvCount_EdgeL     [r] = RecvCount_Fa[r]*3;
         RecvCount_EdgeR     [r] = RecvCount_Fa[r]*3;
         RecvCount_PaddedCr1D[r] = RecvCount_Fa[r];
         RecvCount_MPI_Rank  [r] = RecvCount_Fa[r];
#        ifdef PARTICLE
         RecvCount_NPar      [r] = RecvCount_Fa[r];
#        endif

         RecvDisp_Fa         [r] = ( r == 0 ) ? 0 : RecvDisp_Fa[r-1] + RecvCount_Fa[r-1];
         RecvDisp_Son        [r] = RecvDisp_Fa[r];
         RecvDisp_Sib        [r] = RecvDisp_Fa[r]*26;
         RecvDisp_Cr         [r] = RecvDisp_Fa[r]*3;
         RecvDisp_EdgeL      [r] = RecvDisp_Fa[r]*3;
         RecvDisp_EdgeR      [r] = RecvDisp_Fa[r]*3;
         RecvDisp_PaddedCr1D [r] = RecvDisp_Fa[r];
         RecvDisp_MPI_Rank   [r] = RecvDisp_Fa[r];
#        ifdef PARTICLE
         RecvDisp_NPar       [r] = RecvDisp_Fa[r];
#        endif
      }

//    note that we collect data at one level at a time
      if ( root < 0 ) {
         MPI_Allgatherv( lel.FaList_Local[lv],     amr->NPatchComma[lv][1],      MPI_INT,
                      gel.FaList_AllLv+pc.GID_LvStart[lv],              RecvCount_Fa,           RecvDisp_Fa,         MPI_INT,                   MPI_COMM_WORLD );

         MPI_Allgatherv( lel.SonList_Local[lv],    amr->NPatchComma[lv][1],      MPI_INT,
                      gel.SonList_AllLv+pc.GID_LvStart[lv],             RecvCount_Son,          RecvDisp_Son,        MPI_INT,                   MPI_COMM_WORLD );

         MPI_Allgatherv( lel.SibList_Local[lv][0], amr->NPatchComma[lv][1]*26,   MPI_INT,
                      (gel.SibList_AllLv+pc.GID_LvStart[lv])[0],        RecvCount_Sib,          RecvDisp_Sib,        MPI_INT,                   MPI_COMM_WORLD );

         MPI_Allgatherv( lel.CrList_Local[lv][0],  amr->NPatchComma[lv][1]*3,    MPI_INT,
                      (gel.CrList_AllLv+pc.GID_LvStart[lv])[0],         RecvCount_Cr,           RecvDisp_Cr,         MPI_INT,                   MPI_COMM_WORLD );

         MPI_Allgatherv( lel.EdgeLList_Local[lv][0],  amr->NPatchComma[lv][1]*3, MPI_DOUBLE,
                      (gel.EdgeLList_AllLv+pc.GID_LvStart[lv])[0],      RecvCount_EdgeL,        RecvDisp_EdgeL,      MPI_DOUBLE,                MPI_COMM_WORLD );

         MPI_Allgatherv( lel.EdgeRList_Local[lv][0],  amr->NPatchComma[lv][1]*3, MPI_DOUBLE,
                      (gel.EdgeRList_AllLv+pc.GID_LvStart[lv])[0],      RecvCount_EdgeR,        RecvDisp_EdgeR,      MPI_DOUBLE,                MPI_COMM_WORLD );

         MPI_Allgatherv( lel.PaddedCr1DList_Local[lv],  amr->NPatchComma[lv][1], MPI_UNSIGNED_LONG,
                      (gel.PaddedCr1DList_AllLv+pc.GID_LvStart[lv]),    RecvCount_PaddedCr1D,   RecvDisp_PaddedCr1D, MPI_UNSIGNED_LONG,         MPI_COMM_WORLD );

         MPI_Allgatherv( lel.MPI_RankList_Local[lv],  amr->NPatchComma[lv][1],   MPI_INT,
                      (gel.MPI_RankList_AllLv+pc.GID_LvStart[lv]),      RecvCount_MPI_Rank,     RecvDisp_MPI_Rank,   MPI_INT,                   MPI_COMM_WORLD );

#        ifdef PARTICLE
         MPI_Allgatherv( lel.NParList_Local[lv],    amr->NPatchComma[lv][1],     MPI_INT,
                      gel.NParList_AllLv+pc.GID_LvStart[lv],            RecvCount_NPar,         RecvDisp_NPar,       MPI_INT,                   MPI_COMM_WORLD );
#        endif
      } else {
         MPI_Gatherv( lel.FaList_Local[lv],     amr->NPatchComma[lv][1],         MPI_INT,
                      gel.FaList_AllLv+pc.GID_LvStart[lv],              RecvCount_Fa,           RecvDisp_Fa,         MPI_INT,             root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.SonList_Local[lv],    amr->NPatchComma[lv][1],         MPI_INT,
                      gel.SonList_AllLv+pc.GID_LvStart[lv],             RecvCount_Son,          RecvDisp_Son,        MPI_INT,             root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.SibList_Local[lv][0], amr->NPatchComma[lv][1]*26,      MPI_INT,
                      (gel.SibList_AllLv+pc.GID_LvStart[lv])[0],        RecvCount_Sib,          RecvDisp_Sib,        MPI_INT,             root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.CrList_Local[lv][0],  amr->NPatchComma[lv][1]*3,       MPI_INT,
                      (gel.CrList_AllLv+pc.GID_LvStart[lv])[0],         RecvCount_Cr,           RecvDisp_Cr,         MPI_INT,             root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.EdgeLList_Local[lv][0],  amr->NPatchComma[lv][1]*3,    MPI_DOUBLE,
                      (gel.EdgeLList_AllLv+pc.GID_LvStart[lv])[0],      RecvCount_EdgeL,        RecvDisp_EdgeL,      MPI_DOUBLE,          root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.EdgeRList_Local[lv][0],  amr->NPatchComma[lv][1]*3,    MPI_DOUBLE,
                      (gel.EdgeRList_AllLv+pc.GID_LvStart[lv])[0],      RecvCount_EdgeR,        RecvDisp_EdgeR,      MPI_DOUBLE,          root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.PaddedCr1DList_Local[lv],  amr->NPatchComma[lv][1],    MPI_UNSIGNED_LONG,
                      (gel.PaddedCr1DList_AllLv+pc.GID_LvStart[lv]),    RecvCount_PaddedCr1D,   RecvDisp_PaddedCr1D, MPI_UNSIGNED_LONG,   root, MPI_COMM_WORLD );

         MPI_Gatherv( lel.MPI_RankList_Local[lv],  amr->NPatchComma[lv][1],      MPI_INT,
                      (gel.MPI_RankList_AllLv+pc.GID_LvStart[lv]),      RecvCount_MPI_Rank,     RecvDisp_MPI_Rank,   MPI_INT,             root, MPI_COMM_WORLD );

#        ifdef PARTICLE
         MPI_Gatherv( lel.NParList_Local[lv],    amr->NPatchComma[lv][1],        MPI_INT,
                      gel.NParList_AllLv+pc.GID_LvStart[lv],            RecvCount_NPar,         RecvDisp_NPar,       MPI_INT,             root, MPI_COMM_WORLD );
#        endif
      } // if ( root < 0 ) ... else ...
   } // for (int lv=0; lv<NLEVEL; lv++)

// set isInitialised for ranks that received global exchange lists
   if ( root < 0 ) {
      gel.isInitialised = true;
   } else if ( root == MPI_Rank ) {
      gel.isInitialised = true;
   }

} // FUNCTION : LB_FillGlobalPatchExchangeList



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_ConstructGlobalTree
// Description :  Gather global tree structure as vector indexed by GIDs to root rank
//
// Note        :  - Store global tree AMR structure gathered from all ranks in vector
//                - WARNING: memory allocated for LB_GlobalPatch object must be free by user
//
// Parameter   :  pc   : Reference to LB_PatchCount object
//             :  gel  : Reference to LB_GlobalPatchExchangeList that needs to be initialised by calling LB_FillGlobalPatchExchangeList
//             :  root : Root MPI rank that receives global list, -1 for all ranks
//
// Return      :  - Pointer to LB_GlobalPatch array of length pc.NPatchAllLv allocated on heap
//                - Must be freed by user via delete
//-------------------------------------------------------------------------------------------------------
LB_GlobalPatch* LB_ConstructGlobalTree( LB_PatchCount& pc, LB_GlobalPatchExchangeList& gel, int root ) {

   LB_GlobalPatch* global_tree = NULL;
   if ( root >= 0  &&  root != MPI_Rank )
      return global_tree;

#  ifdef GAMER_DEBUG
   if ( !pc.isInitialised )
      Aux_Error( ERROR_INFO, "call LB_ConstructGlobalTree without initialising LB_PatchCount object !!\n");
   if ( !gel.isInitialised )
      Aux_Error( ERROR_INFO, "call LB_ConstructGlobalTree without initialising LB_GlobalPatchExchangeListobject !!\n");
#  endif

   global_tree = new LB_GlobalPatch[pc.NPatchAllLv];

   int MyGID = 0;
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int i=0; i<NPatchTotal[lv]; ++i) {
         global_tree[MyGID].level      = lv;
         global_tree[MyGID].father     = gel.FaList_AllLv        [MyGID];
         global_tree[MyGID].son        = gel.SonList_AllLv       [MyGID];
         for (int s=0; s<26; s++)
         global_tree[MyGID].sibling[s] = gel.SibList_AllLv       [MyGID][s];
         for (int c=0; c<3 ; c++)
         global_tree[MyGID].corner[c]  = gel.CrList_AllLv        [MyGID][c];
         for (int c=0; c<3 ; c++)
         global_tree[MyGID].EdgeL[c]   = gel.EdgeLList_AllLv     [MyGID][c];
         for (int c=0; c<3 ; c++)
         global_tree[MyGID].EdgeR[c]   = gel.EdgeRList_AllLv     [MyGID][c];
         global_tree[MyGID].PaddedCr1D = gel.PaddedCr1DList_AllLv[MyGID];
         global_tree[MyGID].MPI_Rank   = gel.MPI_RankList_AllLv  [MyGID];
#        ifdef PARTICLE
         global_tree[MyGID].NPar       = gel.NParList_AllLv      [MyGID];
#        endif
         global_tree[MyGID].LB_Idx     = gel.LBIdxList_AllLv     [MyGID];
         MyGID += 1;
      }
   }

   return global_tree;

} // FUNCTION : LB_ConstructGlobalTree



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GatherTree
// Description :  Gather global tree structure as vector indexed by GIDs to root rank
//
// Note        :  - Store global tree AMR structure gathered from all ranks in vector
//                - Initialises pc by calling LB_AllgatherPatchCount ( pc need not be initialised beforehand )
//                - Example usage: Print the information stored in the global tree structure on MPI Node 0
//
//                   LB_PatchCount pc;
//                   LB_GlobalPatch* gt = LB_GatherTree(pc, 0);
//
//                   if ( MPI_Rank == 0 ) {
//                      printf("Information about patches:\n");
//                      for (int i = 0; i < pc.NPatchAllLv; ++i) {
//                         printf("GID %d on level %d residing on MPI rank %d\n", i, gt[i].level, gt[i].MPI_Rank);
//                         printf("Father GID   = %d\n", gt[i].father);
//                         printf("Son GID      = %d\n", gt[i].son);
//                         printf("LB IDx       = %ld\n", gt[i].LB_Idx);
//                         printf("Sibling GIDs = ");
//                         for (int c = 0; c < 26; ++c) {
//                            printf("%d ", gt[i].sibling[c]);
//                         }
//                         printf("\n");
//                         printf("Corners      = ");
//                         for (int c = 0; c < 3; ++c) {
//                            printf("%d ", gt[i].corner[c]);
//                         }
//                         printf("\n");
//                         printf("PaddedCr1D   = %ld\n", gt[i].PaddedCr1D);
//                         printf("Edges:\n");
//                         printf("[x_L, x_R]   = [%9.6f, %9.6f]\n", gt[i].EdgeL[0], gt[i].EdgeR[0]);
//                         printf("[y_L, y_R]   = [%9.6f, %9.6f]\n", gt[i].EdgeL[1], gt[i].EdgeR[1]);
//                         printf("[z_L, z_R]   = [%9.6f, %9.6f]\n", gt[i].EdgeL[2], gt[i].EdgeR[2]);
//
//                         printf("\n");
//#                        ifdef PARTICLE
//                         printf("#Particles   = %d\n", gt[i].NPar);
//#                        endif
//                      }
//                   }
//
//                   delete gt;
//
//                - WARNING: memory allocated for LB_GlobalPatch object must be free by user
//
// Parameter   :  pc   : Reference to LB_PatchCount object
//             :  root : Root MPI rank, -1 for gathering tree at all ranks
//
// Return      :  - Pointer to LB_GlobalPatch array of length pc.NPatchAllLv allocated on heap
//                - Must be freed by user via delete
//-------------------------------------------------------------------------------------------------------
LB_GlobalPatch* LB_GatherTree( LB_PatchCount& pc, int root ) {

// get patch counts per level and per rank from all ranks
   LB_AllgatherPatchCount( pc );

// set up local and global exchange lists
   LB_LocalPatchExchangeList  lel;
   LB_GlobalPatchExchangeList gel( pc, root );

// exchange load balance id's between all ranks
   LB_AllgatherLBIdx( pc, lel, &gel );

// fill local patch lists with information from patches
   LB_FillLocalPatchExchangeList( pc, lel );

// exchange local patch information with other ranks
   LB_FillGlobalPatchExchangeList( pc, lel, gel, root );

// construct and return vector with global tree information
   return LB_ConstructGlobalTree( pc, gel, root );

} // FUNCTION : LB_GatherTree



LB_GlobalTree::LB_GlobalTree(int root) : PatchCount(), Patches(NULL)
{
   Patches = LB_GatherTree(PatchCount, root);
   NPatch  = PatchCount.NPatchAllLv;
}

LB_GlobalTree::~LB_GlobalTree()
{
   delete [] Patches;
}


//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GlobalTree::IsInsidePatch
// Description :  Check whether cell [X, Y, Z] is in patch indexed by GID
//
// Parameter   :  X   : global integer x coordinate, obtained by converting local coordinate with amr->scale
//             :  Y   : global integer y coordinate, obtained by converting local coordinate with amr->scale
//             :  Z   : global integer z coordinate, obtained by converting local coordinate with amr->scale
//             : GID  : global ID of patch
//
// Return      :  true if cell [X, Y, Z] is in patch GID, false otherwise
//-------------------------------------------------------------------------------------------------------
bool LB_GlobalTree::IsInsidePatch(const int X, const int Y, const int Z, const long GID)  const
{
// sanity check
#  ifdef GAMER_DEBUG
   if (GID == -1)
   {
      Aux_Error(ERROR_INFO, "GID == -1!\n");
   }
   if (GID >= NPatch)
   {
      Aux_Error(ERROR_INFO, "GID (%ld) >= NPatch (%ld)!\n", GID, NPatch);
   }
#  endif

   bool IsInside = true;

   int Coordinates[3] = {X, Y, Z};

   for ( int l = 0; l < 3; ++l )
   {
      if (Coordinates[l] < Patches[GID].corner[l] || Coordinates[l] >=  Patches[GID].corner[l] + PS1 * amr->scale[Patches[GID].level])
      {
         IsInside = false;
         break;
      }
   }

   return IsInside;
} // LB_GlobalTree::IsInsidePatch

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GlobalTree::Local2Global
// Description :  Convert local patch coordinates in direction XYZ to global coordinates
//
// Note        :  This function also works for I > PS1 and will return the global coordinate relative to GID.
//                For patch groups, pass the GID of the root patch and the local coordinate I between in {0, ..., PS2}
//                int X = Local2Global(I, 0, GID0);
//                int Y = Local2Global(I, 1, GID0);
//                int Z = Local2Global(I, 2, GID0);
//
// Parameter   :  I   : local coordinate I in {0, ..., PS2} relative to root patch GID
//             : XYZ  : direction
//             : GID  : root patch relative to local coordinate I
//
// Return      :  Global coordinate corresponding to I
//-------------------------------------------------------------------------------------------------------
int LB_GlobalTree::Local2Global(const int I, const int XYZ, const long GID) const
{
// sanity check
#  ifdef GAMER_DEBUG
   if (GID == -1)
   {
      Aux_Error(ERROR_INFO, "GID == -1!\n");
   }
   if (GID >= NPatch)
   {
      Aux_Error(ERROR_INFO, "GID (%ld) >= NPatch (%ld)!\n", GID, NPatch);
   }
#  endif

   return Patches[GID].corner[XYZ] + I * amr->scale[Patches[GID].level];
} // LB_GlobalTree::Local2Global

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GlobalTree::FindRefinedCounterpart
// Description :  Return GID of refined patch with maximum level MaxLv (default = TOP_LEVEL) that global integer coordinates [X, Y, Z] belong to
//
//
// Parameter   :  X   : global integer x coordinate, obtained by converting local coordinate with amr->scale
//             :  Y   : global integer y coordinate, obtained by converting local coordinate with amr->scale
//             :  Z   : global integer z coordinate, obtained by converting local coordinate with amr->scale
//             : GID  : global ID of patch
//             : MaxLv: maximum refinement level
//
// Return      :  GID of patch
//-------------------------------------------------------------------------------------------------------
long LB_GlobalTree::FindRefinedCounterpart(const int X, const int Y, const int Z, const long GID, const int MaxLv)  const
{


#  ifdef GAMER_DEBUG
// sanity check
   if (GID == -1)
   {
      Aux_Error(ERROR_INFO, "GID == -1!\n");
   }
   if (GID >= NPatch)
   {
      Aux_Error(ERROR_INFO, "GID (%ld) >= NPatch (%ld)!\n", GID, NPatch);
   }
#  endif

   if ( Patches[GID].level > MaxLv )
   {
      return -1;
   }


// skip calculation if coordinates are not inside patch GID
   if ( !IsInsidePatch(X, Y, Z, GID) )
   {
      return -1;
   }

   long FaGID  = GID;
   long SonGID = Patches[FaGID].son;

// traverse the tree up until leave nodes
   while ( SonGID != -1 ) {

//    exit loop if FaGID is on MaxLv
      if ( Patches[SonGID].level > MaxLv )
      {
         break;
      }

//    loop over patches in patch group and check if cell {K, J, I} belongs to patch
      for (int LocalID = 0; LocalID < 8; LocalID++)
      {

         if ( IsInsidePatch(X, Y, Z, SonGID + LocalID) )
         {
            FaGID   = SonGID + LocalID;
            SonGID  = Patches[FaGID].son;
            break;

         }
#        ifdef GAMER_DEBUG
         else if ( LocalID == 7 )
         {
            Aux_Error(ERROR_INFO, "Global coordinates {%d, %d, %d} in father patch (GID = %ld, lv = %d), but not in any son patch (GID = %ld, lv = %d)!!\n",
            X, Y, Z, FaGID, Patches[FaGID].level, SonGID, Patches[SonGID].level);
         }
#        endif
      }
   }

#  ifdef GAMER_DEBUG
// sanity check
   if (FaGID == -1)
   {
      Aux_Error(ERROR_INFO, "FaGID == -1!!\n");
   }

// check whether GID is ancestor of FaGID
// number of generations between GID and FaGID
   const int NGenerations = Patches[FaGID].level - Patches[GID].level;

// iterate back through ancestors
   long AncestorGID = FaGID;
   for (int i = 0; i < NGenerations; ++i)
   {
      AncestorGID = Patches[AncestorGID].father;
   }

   if (AncestorGID != GID)
   {
      Aux_Error(ERROR_INFO, "GID (GID = %ld, lv = %d) and Ancestor (%ld) are not related!!\n", GID, Patches[GID].level, AncestorGID);
   }

#  endif

   return FaGID;
} // FUNCTION : LB_GlobalTree::FindRefinedCounterpart


//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GlobalTree::GetPatch
// Description :  Return constant reference to local patch object with GID
//
// Parameter   :  GID  : global ID of patch
//
// Return      :  Constant reference to LB_GlobalPatch
//-------------------------------------------------------------------------------------------------------
const LB_GlobalPatch& LB_GlobalTree::GetPatch(const long GID)  const
{
// sanity check
#  ifdef GAMER_DEBUG
   if (GID == -1)
   {
      Aux_Error(ERROR_INFO, "GID == -1!\n");
   }
   if (GID >= NPatch)
   {
      Aux_Error(ERROR_INFO, "GID (%ld) >= NPatch (%ld)!\n", GID, NPatch);
   }
#  endif


   return Patches[GID];
} // FUNCTION : LB_GlobalTree::GetPatch

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GlobalTree::operator[]
// Description :  Alias for GetPatch
//
// Parameter   :  GID  : global ID of patch
//
// Return      :  Constant reference to LB_GlobalPatch
//-------------------------------------------------------------------------------------------------------
const LB_GlobalPatch& LB_GlobalTree::operator[](const long GID) const
{
   return GetPatch(GID);
} // FUNCTION : LB_GlobalTree::operator[]

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GlobalTree::GetLBPatchCount
// Description :  Return constant reference to LB_PatchCount object
//
// Return      :  Constant reference to LB_PatchCount object
//-------------------------------------------------------------------------------------------------------
const LB_PatchCount&  LB_GlobalTree::GetLBPatchCount() const
{
   return PatchCount;
} // FUNCTION : LB_GlobalTree::GetLBPatchCount
//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GlobalTree::PID2GID
// Description :  Given PID and lv, return GID
//
// Note        :  Does not support patches abroad
//
// Return      :  GID
//-------------------------------------------------------------------------------------------------------
long LB_GlobalTree::PID2GID(const int PID, const int lv) const
{
   if ( PID == -1 )
   {
      return PID;
   }
// patch is a real patch
   else if ( -1 < PID && PID < amr->NPatchComma[lv][1] )
   {
      return PID + PatchCount.GID_Offset[lv];
   }

// patch is a buffer patch (which may lie outside the simulation domain)
   else
   {
#     ifdef GAMER_DEBUG
      Aux_Error( ERROR_INFO, "Lv %d, NPatch %d, PID %d is buffer patch on local MPI Rank %d !!\n",
                 lv, amr->NPatchComma[lv][1], PID, MPI_Rank );
#     endif
      return -1;
   }
} // FUNCTION : LB_GlobalTree::PID2GID
