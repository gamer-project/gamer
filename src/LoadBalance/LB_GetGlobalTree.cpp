#include "GAMER.h"


//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GetGlobalTree
// Description :  Synchronise global tree structure indexed by GIDs to all ranks
//
// Note        :  - Store global tree AMR structure gathered from all ranks in global_tree pointer
//                - Usage: Print level of root patch's son
//                      global_tree_t* t; long t_size; long GID = 0; 
//                      LB_GetGlobalTree(t, t_size);
//                      if (t[GID].son != -1) printf(t[t[GID].son].level);
//                      delete t; 
//                   
//                - User should free global_tree pointer (allocated with new) after usage
//
// Parameter   :  global_tree  : point to newly global tree structure of size NPatchAllLv
//                NPatchAllLv  : length of global_tree array
//                GID_Offset   : array used to store offsets that can be used to convert local PID at level lv to GID via GID = PID + GID_Offset[lv]
//-------------------------------------------------------------------------------------------------------

global_patch_t * LB_GetGlobalTree(long& NPatchAllLv, int GID_Offset[NLEVEL]) {

// pointers to local and global lists 
   long *LBIdxList_Local[NLEVEL], *LBIdxList_AllLv;
   int  (*CrList_Local[NLEVEL])[3], (*CrList_AllLv)[3];
   int  *FaList_Local[NLEVEL], *FaList_AllLv;
   int  *SonList_Local[NLEVEL], *SonList_AllLv;
   int  (*SibList_Local[NLEVEL])[26], (*SibList_AllLv)[26];
#  ifdef PARTICLE
   int  *NParList_Local[NLEVEL], *NParList_AllLv;
#  endif
#  if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
   int  *WaveList_Local[NLEVEL], *WaveList_AllLv;
#  endif // # if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

   long *LBIdxList_Sort[NLEVEL];
   int  *LBIdxList_Sort_IdxTable[NLEVEL];


// temporary variables
   int   MyGID, FaPID, FaGID, FaLv, SonPID, SonGID, SonLv, SibPID, SibGID, MatchIdx;
   long  FaLBIdx, SonLBIdx, SibLBIdx;
   int  *SonCr=NULL, *SibCr=NULL;

// sending and receiving lists for MPI communication 
   int   RecvCount_LBIdx[MPI_NRank], RecvDisp_LBIdx[MPI_NRank], RecvCount_Cr[MPI_NRank], RecvDisp_Cr[MPI_NRank];
   int   RecvCount_Fa[MPI_NRank], RecvDisp_Fa[MPI_NRank], RecvCount_Son[MPI_NRank], RecvDisp_Son[MPI_NRank];
   int   RecvCount_Sib[MPI_NRank], RecvDisp_Sib[MPI_NRank];
#  ifdef PARTICLE
   int   RecvCount_NPar[MPI_NRank], RecvDisp_NPar[MPI_NRank];
#  endif

#  if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
   int   RecvCount_Wave[MPI_NRank], RecvDisp_Wave[MPI_NRank];
#  endif //# if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

// 1. gather the number of patches at different MPI ranks and set the corresponding GID offset
   int (*NPatchAllRank)[NLEVEL] = new int [MPI_NRank][NLEVEL];
   int NPatchLocal[NLEVEL], GID_LvStart[NLEVEL];
   NPatchAllLv=0;

   for (int lv=0; lv<NLEVEL; lv++)  NPatchLocal[lv] = amr->NPatchComma[lv][1];

   MPI_Allgather( NPatchLocal, NLEVEL, MPI_INT, NPatchAllRank[0], NLEVEL, MPI_INT, MPI_COMM_WORLD );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      GID_Offset[lv] = 0;

      for (int r=0; r<MPI_Rank; r++)      GID_Offset[lv] += NPatchAllRank[r][lv];

      for (int FaLv=0; FaLv<lv; FaLv++)   GID_Offset[lv] += NPatchTotal[FaLv];

      NPatchAllLv += NPatchTotal[lv];

      GID_LvStart[lv] = ( lv == 0 ) ? 0 : GID_LvStart[lv-1] + NPatchTotal[lv-1];
   }


// 2. allocate lists 
// global lists for storing tree structure from all ranks
   LBIdxList_AllLv = new long [ NPatchAllLv ];
   CrList_AllLv    = new int  [ NPatchAllLv ][3];
   FaList_AllLv    = new int  [ NPatchAllLv ];
   SonList_AllLv   = new int  [ NPatchAllLv ];
   SibList_AllLv   = new int  [ NPatchAllLv ][26];
#  ifdef PARTICLE
   NParList_AllLv  = new int  [ NPatchAllLv ];
#  endif
#  if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
   WaveList_AllLv  = new int  [ NPatchAllLv ];
#  endif // # if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

// local lists for storing local tree structure
   for (int lv=0; lv<NLEVEL; lv++)
   {
      LBIdxList_Local        [lv] = new long [ amr->NPatchComma[lv][1] ];
      CrList_Local           [lv] = new int  [ amr->NPatchComma[lv][1] ][3];
      FaList_Local           [lv] = new int  [ amr->NPatchComma[lv][1] ];
      SonList_Local          [lv] = new int  [ amr->NPatchComma[lv][1] ];
      SibList_Local          [lv] = new int  [ amr->NPatchComma[lv][1] ][26];
#     ifdef PARTICLE
      NParList_Local         [lv] = new int  [ amr->NPatchComma[lv][1] ];
#     endif
#     if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
      WaveList_Local         [lv] = new int  [ amr->NPatchComma[lv][1] ];
#     endif

      LBIdxList_Sort         [lv] = new long [ NPatchTotal[lv] ];
      LBIdxList_Sort_IdxTable[lv] = new int  [ NPatchTotal[lv] ];
   }


// 3. collect and sort LBIdx from all ranks
// get the search table for converting from PID to GID
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvCount_LBIdx[r] = NPatchAllRank[r][lv];
         RecvDisp_LBIdx [r] = ( r == 0 ) ? 0 : RecvDisp_LBIdx[r-1] + RecvCount_LBIdx[r-1];
      }

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         LBIdxList_Local[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

//    all ranks need to get LBIdxList_Sort since we will use it to calculate GID
      MPI_Allgatherv( LBIdxList_Local[lv], amr->NPatchComma[lv][1], MPI_LONG,
                      LBIdxList_Sort[lv], RecvCount_LBIdx, RecvDisp_LBIdx, MPI_LONG,
                      MPI_COMM_WORLD );
   } // for (int lv=0; lv<NLEVEL; lv++)

// store in the AllLv array BEFORE sorting
   MyGID = 0;

   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<NPatchTotal[lv]; PID++)
      LBIdxList_AllLv[ MyGID++ ] = LBIdxList_Sort[lv][PID];

// sort list and get the corresponding index table (for calculating GID later)
   for (int lv=0; lv<NLEVEL; lv++)
      Mis_Heapsort( NPatchTotal[lv], LBIdxList_Sort[lv], LBIdxList_Sort_IdxTable[lv] );


// 4. store the local tree
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       4-1. LBIdx (set already)
//       LBIdxList_Local[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;


//       4-2. corner
         for (int d=0; d<3; d++)
         CrList_Local[lv][PID][d] = amr->patch[0][lv][PID]->corner[d];


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
            FaGID = FaPID + GID_Offset[FaLv];

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

            Mis_Matching_int( NPatchTotal[FaLv], LBIdxList_Sort[FaLv], 1, &FaLBIdx, &MatchIdx );

#           ifdef GAMER_DEBUG
            if ( MatchIdx < 0 )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d, FaLBIdx %ld, couldn't find a matching patch !!\n",
                       lv, PID, FaPID, FaLBIdx );
#           endif

            FaGID = LBIdxList_Sort_IdxTable[FaLv][MatchIdx] + GID_LvStart[FaLv];
         } // if ( FaPID >= amr->NPatchComma[FaLv][1] )

         FaList_Local[lv][PID] = FaGID;


//       4-4. son GID
         SonPID = amr->patch[0][lv][PID]->son;
         SonLv  = lv + 1;

//       no son (must check this first since SonLv may be out of range --> == NLEVEL)
         if      ( SonPID == -1 )
            SonGID = SonPID;

//       son patch is a real patch at home
         else if ( SonPID >= 0  &&  SonPID < amr->NPatchComma[SonLv][1] )
            SonGID = SonPID + GID_Offset[SonLv];

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

#           if ( defined DEBUG_HDF5  &&  LOAD_BALANCE == HILBERT )
            if ( SonLBIdx - SonLBIdx%8 != 8*amr->patch[0][lv][PID]->LB_Idx )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d, SonCr (%d,%d,%d), incorret SonLBIdx %ld, (MyLBIdx %ld) !!\n",
                       lv, PID, SonPID, SonCr[0], SonCr[1], SonCr[2], SonLBIdx, amr->patch[0][lv][PID]->LB_Idx );
#           endif

            Mis_Matching_int( NPatchTotal[SonLv], LBIdxList_Sort[SonLv], 1, &SonLBIdx, &MatchIdx );

#           ifdef GAMER_DEBUG
            if ( MatchIdx < 0 )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d, SonLBIdx %ld, couldn't find a matching patch !!\n",
                       lv, PID, SonPID, SonLBIdx );
#           endif

            SonGID = LBIdxList_Sort_IdxTable[SonLv][MatchIdx] + GID_LvStart[SonLv];
         } // else if ( SonPID < -1 )

//       son patch is a buffer patch (SonPID >= amr->NPatchComma[SonLv][1]) --> impossible
         else // ( SonPID >= amr->NPatchComma[SonLv][1] )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d is a buffer patch (NRealSonPatch %d) !!\n",
                       lv, PID, SonPID, amr->NPatchComma[SonLv][1] );

         SonList_Local[lv][PID] = SonGID;


//       4-5. sibling GID
         for (int s=0; s<26; s++)
         {
            SibPID = amr->patch[0][lv][PID]->sibling[s];

//          no sibling (SibPID can be either -1 or SIB_OFFSET_NONPERIODIC-BoundaryDirection)
            if      ( SibPID < 0 )
               SibGID = SibPID;

//          sibling patch is a real patch
            else if ( SibPID < amr->NPatchComma[lv][1] )
               SibGID = SibPID + GID_Offset[lv];

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

               Mis_Matching_int( NPatchTotal[lv], LBIdxList_Sort[lv], 1, &SibLBIdx, &MatchIdx );

#              ifdef GAMER_DEBUG
               if ( MatchIdx < 0 )
               Aux_Error( ERROR_INFO, "Lv %d, PID %d, SibPID %d, SibLBIdx %ld, couldn't find a matching patch !!\n",
                          lv, PID, SibPID, SibLBIdx );
#              endif

               SibGID = LBIdxList_Sort_IdxTable[lv][MatchIdx] + GID_LvStart[lv];
            } // if ( SibPID >= amr->NPatchComma[lv][1] )

            SibList_Local[lv][PID][s] = SibGID;

         } // for (int s=0; s<26; s++)


#        ifdef PARTICLE
//       4-6. NPar
         NParList_Local[lv][PID] = amr->patch[0][lv][PID]->NPar;
#        endif

#        if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
//       4-7. Wave flags
         WaveList_Local[lv][PID] = amr->patch[0][lv][PID]->use_wave_flag;
#        endif //#  if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)


// 5. gather data from all ranks
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvCount_Fa  [r] = NPatchAllRank[r][lv];
         RecvCount_Son [r] = RecvCount_Fa[r];
         RecvCount_Sib [r] = RecvCount_Fa[r]*26;
         RecvCount_Cr  [r] = RecvCount_Fa[r]*3;
#        ifdef PARTICLE
         RecvCount_NPar[r] = RecvCount_Fa[r];
#        endif

#        if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
         RecvCount_Wave[r] = RecvCount_Fa[r];
#        endif // #  if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

         RecvDisp_Fa   [r] = ( r == 0 ) ? 0 : RecvDisp_Fa[r-1] + RecvCount_Fa[r-1];
         RecvDisp_Son  [r] = RecvDisp_Fa[r];
         RecvDisp_Sib  [r] = RecvDisp_Fa[r]*26;
         RecvDisp_Cr   [r] = RecvDisp_Fa[r]*3;
#        ifdef PARTICLE
         RecvDisp_NPar [r] = RecvDisp_Fa[r];
#        endif

#        if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
         RecvDisp_Wave[r]  = RecvDisp_Fa[r];
#        endif // #  if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
      }

//    note that we collect data at one level at a time
      MPI_Allgatherv( FaList_Local[lv],     amr->NPatchComma[lv][1],    MPI_INT,
                   FaList_AllLv+GID_LvStart[lv],       RecvCount_Fa,   RecvDisp_Fa,   MPI_INT, MPI_COMM_WORLD );

      MPI_Allgatherv( SonList_Local[lv],    amr->NPatchComma[lv][1],    MPI_INT,
                   SonList_AllLv+GID_LvStart[lv],      RecvCount_Son,  RecvDisp_Son,  MPI_INT, MPI_COMM_WORLD );

      MPI_Allgatherv( SibList_Local[lv][0], amr->NPatchComma[lv][1]*26, MPI_INT,
                   (SibList_AllLv+GID_LvStart[lv])[0], RecvCount_Sib,  RecvDisp_Sib,  MPI_INT, MPI_COMM_WORLD );

      MPI_Allgatherv( CrList_Local[lv][0],  amr->NPatchComma[lv][1]*3,  MPI_INT,
                   (CrList_AllLv+GID_LvStart[lv])[0],  RecvCount_Cr,   RecvDisp_Cr,   MPI_INT, MPI_COMM_WORLD );

#     ifdef PARTICLE
      MPI_Allgatherv( NParList_Local[lv],    amr->NPatchComma[lv][1],    MPI_INT,
                   NParList_AllLv+GID_LvStart[lv],     RecvCount_NPar, RecvDisp_NPar, MPI_INT, MPI_COMM_WORLD );
#     endif

#     if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
      MPI_Allgatherv( WaveList_Local[lv],    amr->NPatchComma[lv][1],    MPI_INT,
                   WaveList_AllLv+GID_LvStart[lv],     RecvCount_Wave, RecvDisp_Wave, MPI_INT, MPI_COMM_WORLD );
#     endif
   } // for (int lv=0; lv<NLEVEL; lv++)


// 6. construct the global tree
   global_patch_t* global_tree = new global_patch_t [NPatchAllLv];

   MyGID = 0;
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int i = 0; i < NPatchTotal[lv]; ++i) {
         global_tree[MyGID].level  = lv; 
         global_tree[MyGID].father = FaList_AllLv[MyGID];
         global_tree[MyGID].son    = SonList_AllLv[MyGID];
         for (int s=0; s<26; s++) global_tree[MyGID].sibling[s] = SibList_AllLv[MyGID][s];
         for (int c=0; c<3 ; c++) global_tree[MyGID].corner[c]  = CrList_AllLv[MyGID][c];
#        ifdef PARTICLE
         global_tree[MyGID].NPar   = NParList_AllLv[MyGID];
#        endif //# ifdef PARTICLE
         global_tree[MyGID].LB_Idx  = LBIdxList_AllLv[MyGID];
         MyGID += 1; 
      }
   }

// 7. free memory

   delete [] NPatchAllRank;

   delete [] LBIdxList_AllLv;
   delete []    CrList_AllLv;
   delete []    FaList_AllLv;
   delete []   SonList_AllLv;
   delete []   SibList_AllLv;
#  ifdef PARTICLE
   delete []  NParList_AllLv;
#  endif

   for (int lv=0; lv<NLEVEL; lv++)
   {
      delete [] LBIdxList_Local[lv];
      delete []    CrList_Local[lv];
      delete []    FaList_Local[lv];
      delete []   SonList_Local[lv];
      delete []   SibList_Local[lv];
#     ifdef PARTICLE
      delete []  NParList_Local[lv];
#     endif

      delete [] LBIdxList_Sort[lv];
      delete [] LBIdxList_Sort_IdxTable[lv];
   }

// return global tree structure
   return global_tree;
}
