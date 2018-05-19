#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_AllocateBufferPatch_Sibling
// Description :  Allocate sibling-buffer patches for all real patches at lv
//
// Note        :  1. This function should be invoked BEFORE LB_AllocateBufferPatch_Father()
//                2. For the base level, an alternative function LB_AllocateBufferPatch_Sibling_Base()
//                   will be invoked (which should be faster)
//                3. During grid refinement, an alternative function LB_Refine_AllocateBufferPatch_Sibling()
//                   will be invoked (which should be faster)
//                4. All buffer patches should be removed in advance
//
// Parameter   :  lv : Target refinement level to allocate sibling-buffer patches
//-------------------------------------------------------------------------------------------------------
void LB_AllocateBufferPatch_Sibling( const int lv )
{

// for the base-level --> invoke the alternative function
   if ( lv == 0 )
   {
      LB_AllocateBufferPatch_Sibling_Base();
      return;
   }


// check the NPatchComma recording
   if ( amr->NPatchComma[lv][1] != amr->num[lv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][1] (%d) != amr->num[%d] (%d)\" !!\n",
                 lv, amr->NPatchComma[lv][1], lv, amr->num[lv] );


   const int PScale        = PATCH_SIZE*amr->scale[lv];
   const int PGScale       = 2*PScale;
   const int NPG_Padded[3] = { amr->BoxScale[0]/PGScale + 2,
                               amr->BoxScale[1]/PGScale + 2,
                               amr->BoxScale[2]/PGScale + 2 };

   int  *Cr0, Cr[3], TRank;
   long  LB_Idx, Coord1D;
   int   MemUnit_Int[MPI_NRank], MemSize_Int[MPI_NRank], MemUnit_Ext[MPI_NRank], MemSize_Ext[MPI_NRank];
   int   MemUnit_Query[MPI_NRank], MemSize_Query[MPI_NRank];
   long *Query_Temp[MPI_NRank], *Int_LBIdx[MPI_NRank], *Ext_Coord1D[MPI_NRank], *Ext_LBIdx_Temp[MPI_NRank];
   long *Ext_LBIdx[MPI_NRank];
   int  *Ext_IdxTable[MPI_NRank];
   int   NQuery_Temp[MPI_NRank], NQuery[MPI_NRank], Int_NQuery_Temp[MPI_NRank], Int_NQuery[MPI_NRank];
   int   Ext_NQuery_Temp[MPI_NRank], Ext_NQuery[MPI_NRank];
   bool  Internal, Allocate;


// 1. prepare send
// ==========================================================================================
// 1.1 set memory size for realloc()
   for (int r=0; r<MPI_NRank; r++)
   {
      MemUnit_Query  [r] = amr->NPatchComma[lv][1];   // set arbitrarily
      MemSize_Query  [r] = MemUnit_Query[r];
      Query_Temp     [r] = (long*)malloc( MemSize_Query[r]*sizeof(long) );
      NQuery_Temp    [r] = 0;
      NQuery         [r] = 0;

      MemUnit_Int    [r] = amr->NPatchComma[lv][1];   // set arbitrarily
      MemSize_Int    [r] = MemUnit_Int[r];
      Int_LBIdx      [r] = (long*)malloc( MemSize_Int[r]*sizeof(long) );
      Int_NQuery_Temp[r] = 0;
      Int_NQuery     [r] = 0;

      MemUnit_Ext    [r] = 100;                       // set arbitrarily
      MemSize_Ext    [r] = MemUnit_Ext[r];
      Ext_Coord1D    [r] = (long*)malloc( MemSize_Ext[r]*sizeof(long) );
      Ext_LBIdx_Temp [r] = (long*)malloc( MemSize_Ext[r]*sizeof(long) );
      Ext_NQuery_Temp[r] = 0;
      Ext_NQuery     [r] = 0;
   }


//###OPTIMIZATION: optimize the treatment of external buffer patches
// 1.2 prepare the LB_Idx query list
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
      Cr0 = amr->patch[0][lv][PID0]->corner;

      for (int k=-1; k<=1; k++)  {  Cr[2] = Cr0[2] + k*PGScale;
      for (int j=-1; j<=1; j++)  {  Cr[1] = Cr0[1] + j*PGScale;
      for (int i=-1; i<=1; i++)  {  Cr[0] = Cr0[0] + i*PGScale;

         LB_Idx = LB_Corner2Index( lv, Cr, CHECK_OFF );

//       ignore LB_Idx lying outside the range
         if ( LB_Idx < amr->LB->CutPoint[lv][0]  ||  LB_Idx >= amr->LB->CutPoint[lv][MPI_NRank] )  continue;

         TRank = LB_Index2Rank( lv, LB_Idx, CHECK_ON );


//       criteria to allocate a buffer patch:
//       --> internal patch: residing in another rank
//           external patch: the external direction is periodic
//       --> internal patches are defined as those with corner[] within the simulation domain
         Internal = true;
         Allocate = true;

         for (int d=0; d<3; d++)
         {
            if ( Cr[d] < 0  ||  Cr[d] >= amr->BoxScale[d] )
            {
               Internal = false;

               if ( OPT__BC_FLU[2*d] != BC_FLU_PERIODIC )
               {
                  Allocate = false;
                  break;
               }
            }
         }

         if ( Internal )   Allocate = ( TRank != MPI_Rank );


//       record the buffer patches to be allocated
         if ( Allocate )
         {
//          a. query list (internal + external buffer patches)
            if ( NQuery_Temp[TRank] >= MemSize_Query[TRank] )
            {
               MemSize_Query[TRank] += MemUnit_Query[TRank];
               Query_Temp   [TRank]  = (long*)realloc( Query_Temp[TRank], MemSize_Query[TRank]*sizeof(long) );
            }

            Query_Temp[TRank][ NQuery_Temp[TRank] ++ ] = LB_Idx;

//          b. internal buffer patch list
            if ( Internal )
            {
               if ( Int_NQuery_Temp[TRank] >= MemSize_Int[TRank] )
               {
                  MemSize_Int[TRank] += MemUnit_Int[TRank];
                  Int_LBIdx  [TRank]  = (long*)realloc( Int_LBIdx[TRank], MemSize_Int[TRank]*sizeof(long) );
               }

               Int_LBIdx[TRank][ Int_NQuery_Temp[TRank] ++ ] = LB_Idx;
            }

//          c. external buffer patch list
            else
            {
               if ( Ext_NQuery_Temp[TRank] >= MemSize_Ext[TRank] )
               {
                  MemSize_Ext   [TRank] += MemUnit_Ext[TRank];
                  Ext_Coord1D   [TRank]  = (long*)realloc( Ext_Coord1D   [TRank],
                                                           MemSize_Ext[TRank]*sizeof(long) );
                  Ext_LBIdx_Temp[TRank]  = (long*)realloc( Ext_LBIdx_Temp[TRank],
                                                           MemSize_Ext[TRank]*sizeof(long) );
               }

               Coord1D = ( (long)(Cr[2]/PGScale+1)*NPG_Padded[1] + Cr[1]/PGScale+1 )*NPG_Padded[0]
                         + Cr[0]/PGScale+1;
               Ext_Coord1D    [TRank][ Ext_NQuery_Temp[TRank] ] = Coord1D;
               Ext_LBIdx_Temp [TRank][ Ext_NQuery_Temp[TRank] ] = LB_Idx;
               Ext_NQuery_Temp[TRank] ++;
            }

         } // if ( Allocate )
      }}} // k,j,i
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


// 1.3 sort the query list and remove duplicates (with the same LB_Idx)
   for (int r=0; r<MPI_NRank; r++)
   {
      Mis_Heapsort( NQuery_Temp[r], Query_Temp[r], NULL );

      if ( NQuery_Temp[r] > 0 )  NQuery[r] = 1;

      for (int t=1; t<NQuery_Temp[r]; t++)
         if ( Query_Temp[r][t] != Query_Temp[r][t-1] )   Query_Temp[r][ NQuery[r]++ ] = Query_Temp[r][t];
   }


// 1.4 sort the internal buffer patch list and remove duplicates (with the same LB_Idx)
   for (int r=0; r<MPI_NRank; r++)
   {
      Mis_Heapsort( Int_NQuery_Temp[r], Int_LBIdx[r], NULL );

      if ( Int_NQuery_Temp[r] > 0 )  Int_NQuery[r] = 1;

      for (int t=1; t<Int_NQuery_Temp[r]; t++)
         if ( Int_LBIdx[r][t] != Int_LBIdx[r][t-1] )   Int_LBIdx[r][ Int_NQuery[r] ++ ] = Int_LBIdx[r][t];
   }


// 1.5 sort the external buffer patch list and remove duplicates (with the same 1D coord)
// ***note that different external buffer patches can map to the same internal real patches in periodic B.C.
   for (int r=0; r<MPI_NRank; r++)
   {
      Ext_IdxTable[r] = (int*)malloc( Ext_NQuery_Temp[r]*sizeof(int) );

      Mis_Heapsort( Ext_NQuery_Temp[r], Ext_Coord1D[r], Ext_IdxTable[r] );

      if ( Ext_NQuery_Temp[r] > 0 )  Ext_NQuery[r] = 1;

      for (int t=1; t<Ext_NQuery_Temp[r]; t++)
      {
         if ( Ext_Coord1D[r][t] != Ext_Coord1D[r][t-1] )
         {
            Ext_Coord1D [r][ Ext_NQuery[r] ] = Ext_Coord1D [r][t];
            Ext_IdxTable[r][ Ext_NQuery[r] ] = Ext_IdxTable[r][t];
            Ext_NQuery[r] ++;
         }
      }

//    construct the final mapping list for the external buffer patches
      Ext_LBIdx[r] = (long*)malloc( Ext_NQuery[r]*sizeof(long) );

      for (int t=0; t<Ext_NQuery[r]; t++)    Ext_LBIdx[r][t] = Ext_LBIdx_Temp[r][ Ext_IdxTable[r][t] ];

      Mis_Heapsort( Ext_NQuery[r], Ext_LBIdx[r], Ext_IdxTable[r] );
   }



// 2. transfer data : (SendBuf_Query --> RecvBuf_Query --> SendBuf_Reply --> RecvBuf_Reply)
// ==========================================================================================
   int   Query_Disp[MPI_NRank], Reply_Disp[MPI_NRank], NReply[MPI_NRank], NQuery_Total, NReply_Total, Counter;
   long *SendBuf_Query, *RecvBuf_Query;
   char *SendBuf_Reply, *RecvBuf_Reply;

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


// 2.4 prepare answers
   for (int r=0; r<MPI_NRank; r++)
      Mis_Matching_char( amr->NPatchComma[lv][1], amr->LB->IdxList_Real[lv], NReply[r],
                         RecvBuf_Query+Reply_Disp[r], SendBuf_Reply+Reply_Disp[r] );


// 2.5 send replies
   MPI_Alltoallv( SendBuf_Reply, NReply, Reply_Disp, MPI_CHAR,
                  RecvBuf_Reply, NQuery, Query_Disp, MPI_CHAR, MPI_COMM_WORLD );



// 3. allocate buffer patches
// ==========================================================================================
   int  Int_Counter, Ext_Counter, Counter2=0;
   long Int_Target, Ext_Target, Cr1D;

   for (int r=0; r<MPI_NRank; r++)
   {
      Int_Counter = 0;
      Ext_Counter = 0;
      Int_Target  = ( Int_Counter < Int_NQuery[r] ) ? Int_LBIdx[r][Int_Counter] : -1;
      Ext_Target  = ( Ext_Counter < Ext_NQuery[r] ) ? Ext_LBIdx[r][Ext_Counter] : -1;

      for (int t=0; t<NQuery[r]; t++)
      {
         LB_Idx = Query_Temp[r][t];

//       validate LB_Idx
#        ifdef GAMER_DEBUG
         if ( LB_Idx != Int_Target  &&  LB_Idx != Ext_Target )
            Aux_Error( ERROR_INFO, "incorrect LB_Idx (%ld) != Int_Target (%ld) and Ext_Target (%ld) !!\n",
                       LB_Idx, Int_Target, Ext_Target );
#        endif

         if ( RecvBuf_Reply[ Counter2 ++ ] == 1 )
         {
//          note that a given LB_Idx can map to both internal and external patches at the same time
            if ( LB_Idx == Int_Target ) // internal buffer patches
            {
               LB_Index2Corner( lv, LB_Idx, Cr, CHECK_ON );

//             father patch is still unkown, data array is not allocated yet
               amr->pnew( lv, Cr[0],        Cr[1],        Cr[2],        -1, false, false, false );
               amr->pnew( lv, Cr[0]+PScale, Cr[1],        Cr[2],        -1, false, false, false );
               amr->pnew( lv, Cr[0],        Cr[1]+PScale, Cr[2],        -1, false, false, false );
               amr->pnew( lv, Cr[0],        Cr[1],        Cr[2]+PScale, -1, false, false, false );
               amr->pnew( lv, Cr[0]+PScale, Cr[1]+PScale, Cr[2],        -1, false, false, false );
               amr->pnew( lv, Cr[0],        Cr[1]+PScale, Cr[2]+PScale, -1, false, false, false );
               amr->pnew( lv, Cr[0]+PScale, Cr[1],        Cr[2]+PScale, -1, false, false, false );
               amr->pnew( lv, Cr[0]+PScale, Cr[1]+PScale, Cr[2]+PScale, -1, false, false, false );

               amr->NPatchComma[lv][2] += 8;

               Int_Counter ++;
               Int_Target = ( Int_Counter < Int_NQuery[r] ) ? Int_LBIdx[r][Int_Counter] : -1;
            }

            if ( LB_Idx == Ext_Target ) // external buffer patches
            {
//             loop over all external buffer patches mapping to the same internal real patches
               do
               {
                  Cr1D  = Ext_Coord1D[r][ Ext_IdxTable[r][Ext_Counter] ];
                  Cr[0] = PGScale*( -1 + Cr1D%NPG_Padded[0]                               );
                  Cr[1] = PGScale*( -1 + Cr1D%(NPG_Padded[0]*NPG_Padded[1])/NPG_Padded[0] );
                  Cr[2] = PGScale*( -1 + Cr1D/(NPG_Padded[0]*NPG_Padded[1])               );

//                father patch is still unkown, data array is not allocated yet
                  amr->pnew( lv, Cr[0],        Cr[1],        Cr[2],        -1, false, false, false );
                  amr->pnew( lv, Cr[0]+PScale, Cr[1],        Cr[2],        -1, false, false, false );
                  amr->pnew( lv, Cr[0],        Cr[1]+PScale, Cr[2],        -1, false, false, false );
                  amr->pnew( lv, Cr[0],        Cr[1],        Cr[2]+PScale, -1, false, false, false );
                  amr->pnew( lv, Cr[0]+PScale, Cr[1]+PScale, Cr[2],        -1, false, false, false );
                  amr->pnew( lv, Cr[0],        Cr[1]+PScale, Cr[2]+PScale, -1, false, false, false );
                  amr->pnew( lv, Cr[0]+PScale, Cr[1],        Cr[2]+PScale, -1, false, false, false );
                  amr->pnew( lv, Cr[0]+PScale, Cr[1]+PScale, Cr[2]+PScale, -1, false, false, false );

                  amr->NPatchComma[lv][2] += 8;

                  Ext_Counter ++;
                  Ext_Target = ( Ext_Counter < Ext_NQuery[r] ) ? Ext_LBIdx[r][Ext_Counter] : -1;

//                check : no external buffer patches should be allocated for the non-periodic B.C.
#                 ifdef GAMER_DEBUG
                  for (int d=0; d<3; d++)
                  {
                     if (  OPT__BC_FLU[2*d] != BC_FLU_PERIODIC  &&  ( Cr[d] < 0 || Cr[d] >= amr->BoxScale[d] )  )
                        Aux_Error( ERROR_INFO, "Cr[%d] = %d lies outside the simulation box for non-periodic BC !!\n",
                                   d, Cr[d] );
                  }
#                 endif
               }
               while ( Ext_Target == Ext_LBIdx[r][Ext_Counter-1] );

            } // if ( LB_Idx != Ext_LBIdx[r][Ext_Counter] ) ... else ...
         } // if ( RecvBuf_Reply[ Counter2 ++ ] == 1 )

         else
         {
//          note that a given LB_Idx can map to both internal and external patches at the same time
            if ( LB_Idx == Int_Target ) // skip the non-existing internal buffer patches
            {
               Int_Counter ++;
               Int_Target = ( Int_Counter < Int_NQuery[r] ) ? Int_LBIdx[r][Int_Counter] : -1;
            }

            if ( LB_Idx == Ext_Target ) // skip the non-existing external buffer patches
            {
               do
               {
                  Ext_Counter ++;
                  Ext_Target = ( Ext_Counter < Ext_NQuery[r] ) ? Ext_LBIdx[r][Ext_Counter] : -1;
               }
               while ( Ext_Target == Ext_LBIdx[r][Ext_Counter-1] );
            }
         } // if ( RecvBuf_Reply[ Counter2++ ] == 1 ) ... else ...

      } // for (int t=0; t<NQuery[r]; t++)
   } // for (int r=0; r<MPI_NRank; r++)


// record NPatchComma
   for (int m=3; m<28; m++)   amr->NPatchComma[lv][m] = amr->NPatchComma[lv][2];

   if ( amr->NPatchComma[lv][2] != amr->num[lv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][2] (%d) != amr->num[%d] (%d)\" !!\n",
                 lv, amr->NPatchComma[lv][2], lv, amr->num[lv] );



// 4. record the padded 1D corner coordinates (which can be overwritten by "LB_AllocateBufferPatch_Father")
// ==========================================================================================
   const int NPatch = amr->NPatchComma[lv][2];

   amr->LB->PaddedCr1DList         [lv] = (ulong*)realloc( amr->LB->PaddedCr1DList         [lv],
                                                           NPatch*sizeof(ulong) );
   amr->LB->PaddedCr1DList_IdxTable[lv] = (int*  )realloc( amr->LB->PaddedCr1DList_IdxTable[lv],
                                                           NPatch*sizeof(int ) );

   for (int PID=0; PID<NPatch; PID++)
      amr->LB->PaddedCr1DList[lv][PID] = amr->patch[0][lv][PID]->PaddedCr1D;

   Mis_Heapsort( NPatch, amr->LB->PaddedCr1DList[lv], amr->LB->PaddedCr1DList_IdxTable[lv] );

// check : no duplicate patches at lv
#  ifdef GAMER_DEBUG
   for (int t=1; t<amr->num[lv]; t++)
   {
      if ( amr->LB->PaddedCr1DList[lv][t] == amr->LB->PaddedCr1DList[lv][t-1] )
         Aux_Error( ERROR_INFO, "duplicate patches at lv %d, PaddedCr1D %lu, PID = %d and %d !!\n",
                    lv, amr->LB->PaddedCr1DList[lv][t], amr->LB->PaddedCr1DList_IdxTable[lv][t],
                    amr->LB->PaddedCr1DList_IdxTable[lv][t-1] );
   }
#  endif


// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( Query_Temp    [r] );
      free( Int_LBIdx     [r] );
      free( Ext_Coord1D   [r] );
      free( Ext_IdxTable  [r] );
      free( Ext_LBIdx     [r] );
      free( Ext_LBIdx_Temp[r] );
   }
   delete [] SendBuf_Query;
   delete [] RecvBuf_Query;
   delete [] SendBuf_Reply;
   delete [] RecvBuf_Reply;

} // FUNCTION : LB_AllocateBufferPatch_Sibling



#endif // #ifdef LOAD_BALANCE
