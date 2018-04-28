#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_FindSonNotHome
// Description :  Assign son indices for father patches without sons at home (these father patches can either
//                have no sons or have sons not home)
//
// Note        :  1. Apply to both real and buffer patches
//                2. For patches with sons living abroad (not at the same rank as iteself), their son indices
//                   are set to "SON_OFFSET_LB-SonRank", where SonRank represents the MPI rank where their sons live.
//                   (e.g., if the real son is at rank 123, the son index is set to SON_OFFSET_LB-123=-1123 for
//                   SON_OFFSET_LB==-1000)
//                   --> All patches with sons will have SonPID != -1
//                3. This function will search over all father patches with SonPID <= -1
//                4. Son indices will always be real patches at SonLv (or SON_OFFSET_LB-SonRank, where SonRank
//                   is the home rank of the **real** son patch)
//                   --> Even if there is an external buffer patch at FaLv with the same PaddedCr1D as an
//                       external buffer patch at SonLv and both the two buffer patches and the real patch
//                       corresponding to the buffer patch at SonLv are in the same rank, the son index of the
//                       buffer patch at FaLv will still be set to "SON_OFFSET_LB-SonRank" with SonRank == MPI_Rank
//
// Parameter   :  FaLv        : Target refinement level of fathers
//                SearchAllFa : Whether to search over all father patches or not
//                NInput      : Number of target father patches in "TargetFaPID"
//                              (useful only if "SearchAllFa == false")
//                TargetFaPID : Lists recording all target father patches
//                              (useful only if "SearchAllFa == false")
//-------------------------------------------------------------------------------------------------------
void LB_FindSonNotHome( const int FaLv, const bool SearchAllFa, const int NInput, int* TargetFaPID )
{

// nothing to do for the maximum level
   if ( FaLv == NLEVEL - 1 )  return;


// check
   if ( FaLv < 0  ||  FaLv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaLv", FaLv );

   if ( !SearchAllFa  &&  NInput != 0  &&  TargetFaPID == NULL )
      Aux_Error( ERROR_INFO, "FaLv %d, NInput %d, TargetFaPID == NULL !!\n", FaLv, NInput );


   const int  SonLv     = FaLv + 1;
   const int  NFaPatch  = amr->num[FaLv];    // real + buffer patches
   const int  NTargetFa = ( SearchAllFa ) ? NFaPatch : NInput;

   int FaPID;


// 0. construct the target father patch list
// ==========================================================================================
   if ( SearchAllFa )
   {
      TargetFaPID = new int [NTargetFa];

      for (int t=0; t<NTargetFa; t++)  TargetFaPID[t] = t;
   }


// 1 construct the query list for different ranks
// ==========================================================================================
   int   TRank, MemUnit_Query[MPI_NRank], MemSize_Query[MPI_NRank], NQuery[MPI_NRank];
   int  *Cr, *FaPID_List[MPI_NRank], *FaPID_IdxTable[MPI_NRank];
   long  LB_Idx;
   long *Query_Temp[MPI_NRank];
   bool  Internal;

// 1.1 set memory allocation unit
   for (int r=0; r<MPI_NRank; r++)
   {
      MemUnit_Query[r] = 1 + NTargetFa/MPI_NRank;  // set arbitrarily
      MemSize_Query[r] = MemUnit_Query[r];
      Query_Temp   [r] = (long*)malloc( MemSize_Query[r]*sizeof(long) );
      FaPID_List   [r] = (int* )malloc( MemSize_Query[r]*sizeof(int ) );
      NQuery       [r] = 0;
   }

// 1.2 construct the unsorted query list for different ranks
   for (int t=0; t<NTargetFa; t++)
   {
      FaPID = TargetFaPID[t];

      if ( amr->patch[0][FaLv][FaPID]->son <= -1 )
      {
         Cr = amr->patch[0][FaLv][FaPID]->corner;

//###NOTE: faster version can only be applied to the Hilbert space-filling curve
#        if ( LOAD_BALANCE == HILBERT )
         LB_Idx = 8*amr->patch[0][FaLv][FaPID]->LB_Idx;     // faster, LB_Idx of one of the eight sons
#        else
         LB_Idx = LB_Corner2Index( SonLv, Cr, CHECK_OFF );  // LB_Idx of son 0
#        endif

//       ignore LB_Idx lying outside the range
         if ( LB_Idx < amr->LB->CutPoint[SonLv][0]  ||  LB_Idx >= amr->LB->CutPoint[SonLv][MPI_NRank] )
         {
            amr->patch[0][FaLv][FaPID]->son = -1;
            continue;
         }

//       ignore internal patches in the same rank
         TRank    = LB_Index2Rank( SonLv, LB_Idx, CHECK_ON );
         Internal = true;

         for (int d=0; d<3; d++)
         {
            if ( Cr[d] < 0  ||  Cr[d] >= amr->BoxScale[d] )
            {
               Internal = false;

#              ifdef GAMER_DEBUG
               if ( OPT__BC_FLU[2*d] != BC_FLU_PERIODIC )
                  Aux_Error( ERROR_INFO, "external patch (lv %d, PID %d) for the non-periodic BC !!\n", FaLv, FaPID );
#              endif

               break;
            }
         }

         if ( Internal  &&  TRank == MPI_Rank )
         {
#           ifdef GAMER_DEBUG
            if ( amr->patch[0][FaLv][FaPID]->son < -1 )
               Aux_Error( ERROR_INFO, "FaLv %d, FaPID (%d)'s SonPID < -1 !!\n",
                          FaLv, FaPID, amr->patch[0][FaLv][FaPID]->son );
#           endif

            continue;
         }

//       allocate enough memory
         if ( NQuery[TRank] >= MemSize_Query[TRank] )
         {
            MemSize_Query[TRank] += MemUnit_Query[TRank];
            Query_Temp   [TRank]  = (long*)realloc( Query_Temp[TRank], MemSize_Query[TRank]*sizeof(long) );
            FaPID_List   [TRank]  = (int* )realloc( FaPID_List[TRank], MemSize_Query[TRank]*sizeof(int ) );
         }

//       record list
         Query_Temp[TRank][ NQuery[TRank] ] = LB_Idx;
         FaPID_List[TRank][ NQuery[TRank] ] = FaPID;
         NQuery    [TRank] ++;

      } // if ( amr->patch[0][FaLv][FaPID]->son <= -1 )
   } // for (int t=0; t<NTargetFa; t++)

// 1.3 sort the query list for different ranks
   for (int r=0; r<MPI_NRank; r++)
   {
      FaPID_IdxTable[r] = new int [ NQuery[r] ];

      Mis_Heapsort( NQuery[r], Query_Temp[r], FaPID_IdxTable[r] );
   }


// 2 transfer data : (SendBuf_Query --> RecvBuf_Query --> SendBuf_Reply --> RecvBuf_Reply)
// ==========================================================================================
   int   Query_Disp[MPI_NRank], Reply_Disp[MPI_NRank], NReply[MPI_NRank], NQuery_Total, NReply_Total, Counter;
   long *SendBuf_Query=NULL, *RecvBuf_Query=NULL;
   char *SendBuf_Reply=NULL, *RecvBuf_Reply=NULL;

// 2.1 send the number of queries
   MPI_Alltoall( NQuery, 1, MPI_INT, NReply, 1, MPI_INT, MPI_COMM_WORLD );

// 2.2 prepare the query and reply arrays
   Query_Disp[0] = 0;
   Reply_Disp[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Query_Disp[r] = Query_Disp[r-1] + NQuery[r-1];
      Reply_Disp[r] = Reply_Disp[r-1] + NReply[r-1];
   }
   NQuery_Total = Query_Disp[MPI_NRank-1] + NQuery[MPI_NRank-1];
   NReply_Total = Reply_Disp[MPI_NRank-1] + NReply[MPI_NRank-1];

   SendBuf_Query = new long [NQuery_Total];
   RecvBuf_Query = new long [NReply_Total];
   SendBuf_Reply = new char [NReply_Total];
   RecvBuf_Reply = new char [NQuery_Total];

   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<NQuery[r]; t++)
      SendBuf_Query[ Counter ++ ] = Query_Temp[r][t];

// 2.3 send queries
   MPI_Alltoallv( SendBuf_Query, NQuery, Query_Disp, MPI_LONG,
                  RecvBuf_Query, NReply, Reply_Disp, MPI_LONG, MPI_COMM_WORLD );

// 2.4 prepare replies
   for (int r=0; r<MPI_NRank; r++)
      Mis_Matching_char( amr->NPatchComma[SonLv][1], amr->LB->IdxList_Real[SonLv], NReply[r],
                         RecvBuf_Query+Reply_Disp[r], SendBuf_Reply+Reply_Disp[r] );

// 2.5 send replies
   MPI_Alltoallv( SendBuf_Reply, NReply, Reply_Disp, MPI_CHAR,
                  RecvBuf_Reply, NQuery, Query_Disp, MPI_CHAR, MPI_COMM_WORLD );


// 3 set SonPID to "SON_OFFSET_LB-SonRank"
// ==========================================================================================
   Counter=0;

   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<NQuery[r]; t++)
   {
      FaPID = FaPID_List[r][ FaPID_IdxTable[r][t] ];

      if ( RecvBuf_Reply[ Counter ++ ] == 1 )
      {
         amr->patch[0][FaLv][FaPID]->son = SON_OFFSET_LB - r;

//       check : only external buffer patches can have sons at home but with son indices < -1
#        ifdef GAMER_DEBUG
         Cr = amr->patch[0][FaLv][FaPID]->corner;

         Internal = true;
         for (int d=0; d<3; d++)
         {
            if ( Cr[d] < 0  ||  Cr[d] >= amr->BoxScale[d] )
            {
               Internal = false;
               break;
            }
         }

         if ( Internal  &&  r == MPI_Rank )
            Aux_Error( ERROR_INFO, "FaLv %d, FaPID %d's son should be home !!\n", FaLv, FaPID );
#        endif
      } // if ( RecvBuf_Reply[ Counter ++ ] == 1 )

      else
         amr->patch[0][FaLv][FaPID]->son = -1;

   } // for (int r=0; r<MPI_NRank; r++) ... for (int t=0; t<NQuery[r]; t++)


// 4. check results in debug mode
#  ifdef GAMER_DEBUG
   const int SonNReal = amr->NPatchComma[SonLv][1];
   int SonPID, SonPID0;

// check 1 : all real patches must have fathers
   for (SonPID=0; SonPID<SonNReal; SonPID++)
   {
      if ( amr->patch[0][SonLv][SonPID]->father == -1 )
         Aux_Error( ERROR_INFO, "Check 1, SonLv %d: SonPID (%d) has no father !!\n", SonLv, SonPID );
   }

// check 2 : SonPIDs should always be real patches at SonLv
   for (FaPID=0; FaPID<NFaPatch; FaPID++)
   {
      SonPID = amr->patch[0][FaLv][FaPID]->son;

      if ( SonPID >= SonNReal )
         Aux_Error( ERROR_INFO, "Check 2, FaLv %d: FaPID (%d) -> son (%d) >= SonNReal (%d) !!\n",
                    FaLv, FaPID, SonPID, SonNReal );
   }

// check 3 : father's son = itself
   for (SonPID=0; SonPID<SonNReal; SonPID++)
   {
      FaPID   = amr->patch[0][SonLv][SonPID]->father;
      SonPID0 = SonPID - SonPID%8;

      if ( amr->patch[0][FaLv][FaPID]->son != SonPID0 )
         Aux_Error( ERROR_INFO, "Check 3, SonLv %d: FaPID (%d) -> son (%d) != SonPID (%d) -> SonPID0 (%d) !!\n",
                    SonLv, FaPID, amr->patch[0][FaLv][FaPID]->son, SonPID, SonPID0 );
   }

// check 4 : son's father = itself
   for (FaPID=0; FaPID<NFaPatch; FaPID++)
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

   for (SonPID0=0; SonPID0<SonNReal; SonPID0+=8)
   {
      Corner_Son = amr->patch[0][SonLv][SonPID0]->corner;

      for (FaPID=0; FaPID<NFaPatch; FaPID++)
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
#  endif // #ifdef GAMER_DEBUG


// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( Query_Temp[r] );
      free( FaPID_List[r] );
      delete [] FaPID_IdxTable[r];
   }

   delete [] SendBuf_Query;
   delete [] RecvBuf_Query;
   delete [] SendBuf_Reply;
   delete [] RecvBuf_Reply;
   if ( SearchAllFa )   delete [] TargetFaPID;

} // FUNCTION : LB_FindSonNotHome



#endif // #ifdef LOAD_BALANCE
