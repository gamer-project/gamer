#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_LB_Refine_GetCoarseFineInterfaceBField
// Description :  Get the magnetic field on the coarse-fine interfaces for the divergence-free interpolation
//
// Note        :  1. Invoked by LB_Refine()
//                2. Must invoke LB_Refine_GetNewRealPatchList() first to prepare the requisite info
//                3. Divergence-free interpolation is conducted by LB_Refine_AllocateNewPatch()
//
// Parameter   :  FaLv              : Parent level
//                NNew_Home         : Number of home patches at FaLv to allocate son patches
//                NNew_Away         : Number of away patches at FaLv to allocate son patches
//                CFB_SibLBIdx_Home : Load-balance indices of the siblings of home patches
//                CFB_SibLBIdx_Away : Load-balance indices of the siblings of away patches
//                CFB_SibRank_Home  : MPI ranks of the siblings of home patches
//                CFB_SibRank_Away  : MPI ranks of the siblings of away patches
//                CFB_BField        : Coarse-fine interface B field array
//                CFB_NSibEachRank  : Number of siblings receiving data from each rank (including its own rank)
//
// Return      :  CFB_SibRank_Home, CFB_SibRank_Away, CFB_BField, CFB_NSibEachRank
//-------------------------------------------------------------------------------------------------------
void MHD_LB_Refine_GetCoarseFineInterfaceBField(
   const int FaLv, const int NNew_Home, const int NNew_Away,
   const long (*CFB_SibLBIdx_Home)[6], const long (*CFB_SibLBIdx_Away)[6],
   int (*&CFB_SibRank_Home)[6], int (*&CFB_SibRank_Away)[6],
   real *&CFB_BField, long *CFB_NSibEachRank )
{

   const int SonLv   = FaLv + 1;
   const int MemUnit = 1 + ( NNew_Home + NNew_Away )/MPI_NRank;

   int   MemSize              [MPI_NRank];
   long  RecvEachRank_N       [MPI_NRank];
   int  *RecvEachRank_SibID   [MPI_NRank];
   long *RecvEachRank_SibLBIdx[MPI_NRank];
   long  SendEachRank_N       [MPI_NRank];
   int  *SendEachRank_SibID   [MPI_NRank];
   int  *SendEachRank_PID     [MPI_NRank];

   for (int r=0; r<MPI_NRank; r++)
   {
      MemSize              [r] = MemUnit;
      RecvEachRank_N       [r] = 0L;
      RecvEachRank_SibID   [r] = (int* )malloc( MemSize[r]*sizeof(int ) );
      RecvEachRank_SibLBIdx[r] = (long*)malloc( MemSize[r]*sizeof(long) );
   }

   CFB_SibRank_Home = new int [NNew_Home][6];
   CFB_SibRank_Away = new int [NNew_Away][6];


// 1. set the parameters for receiving fine-grid B field from other ranks
//    --> need to check both home and away patches
   const int  NNew_Both[2]               = { NNew_Home, NNew_Away };
   const long (*CFB_SibLBIdx_Both[2])[6] = { CFB_SibLBIdx_Home, CFB_SibLBIdx_Away };
         int  (*CFB_SibRank_Both[2])[6]  = { CFB_SibRank_Home, CFB_SibRank_Away };

   for (int g=0; g<2; g++)
   for (int t=0; t<NNew_Both[g]; t++)
   for (int s=0; s<6; s++)
   {
      const long SibLBIdx = CFB_SibLBIdx_Both[g][t][s];

//    SibLBIdx is set to -1 in LB_Refine_GetNewRealPatchList.cpp() on the coarse-coarse interfaces
      if ( SibLBIdx >= 0 )
      {
         const int TRank = LB_Index2Rank( SonLv, SibLBIdx, CHECK_ON );

         CFB_SibRank_Both[g][t][s] = TRank;

         if ( RecvEachRank_N[TRank] >= (long)MemSize[TRank] )
         {
             MemSize              [TRank] += MemUnit;
             RecvEachRank_SibID   [TRank]  = (int* )realloc( RecvEachRank_SibID   [TRank], MemSize[TRank]*sizeof(int ) );
             RecvEachRank_SibLBIdx[TRank]  = (long*)realloc( RecvEachRank_SibLBIdx[TRank], MemSize[TRank]*sizeof(long) );
         }

         RecvEachRank_SibID   [TRank][ RecvEachRank_N[TRank] ] = s;
         RecvEachRank_SibLBIdx[TRank][ RecvEachRank_N[TRank] ] = SibLBIdx;
         RecvEachRank_N       [TRank] ++;
      } // if ( SibLBIdx >= 0 )

//    set TRank to a negative value to indicate coarse-coarse interfaces
      else
         CFB_SibRank_Both[g][t][s] = -1;
   } // for g, t, s



// 2. set the parameters for sending fine-grid B field to other ranks
   int   Counter;
   long  Send_NTotal, Recv_NTotal;
   long  Send_NList[MPI_NRank], Recv_NList[MPI_NRank];
   long  Send_Disp[MPI_NRank], Recv_Disp[MPI_NRank];
   int  *SendBuf_SibID=NULL, *RecvBuf_SibID=NULL;
   long *SendBuf_LBIdx=NULL, *RecvBuf_LBIdx=NULL;

// 2.1. prepare the send info
   MPI_Alltoall( RecvEachRank_N, 1, MPI_LONG, SendEachRank_N, 1, MPI_LONG, MPI_COMM_WORLD );

   memcpy( Send_NList, RecvEachRank_N, MPI_NRank*sizeof(long) );
   memcpy( Recv_NList, SendEachRank_N, MPI_NRank*sizeof(long) );

   Send_Disp[0] = 0L;
   Recv_Disp[0] = 0L;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp[r] = Send_Disp[r-1] + Send_NList[r-1];
      Recv_Disp[r] = Recv_Disp[r-1] + Recv_NList[r-1];
   }
   Send_NTotal = Send_Disp[MPI_NRank-1] + Send_NList[MPI_NRank-1];
   Recv_NTotal = Recv_Disp[MPI_NRank-1] + Recv_NList[MPI_NRank-1];

   SendBuf_SibID = new int  [Send_NTotal];
   RecvBuf_SibID = new int  [Recv_NTotal];
   SendBuf_LBIdx = new long [Send_NTotal];
   RecvBuf_LBIdx = new long [Recv_NTotal];

   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<Send_NList[r]; t++)
   {
      SendBuf_SibID[Counter] = RecvEachRank_SibID   [r][t];
      SendBuf_LBIdx[Counter] = RecvEachRank_SibLBIdx[r][t];
      Counter ++;
   }


// 2.2. send --> recv
   MPI_Alltoallv_GAMER( SendBuf_SibID, Send_NList, Send_Disp, MPI_INT,
                        RecvBuf_SibID, Recv_NList, Recv_Disp, MPI_INT,  MPI_COMM_WORLD );

   MPI_Alltoallv_GAMER( SendBuf_LBIdx, Send_NList, Send_Disp, MPI_LONG,
                        RecvBuf_LBIdx, Recv_NList, Recv_Disp, MPI_LONG, MPI_COMM_WORLD );


// 2.3. record the received info
   for (int r=0; r<MPI_NRank; r++)
   {
//    2.3.1. find the matching patches for sending fine-grid B field
      const int  *RecvPtr_SibID = RecvBuf_SibID + Recv_Disp[r];
            long *RecvPtr_LBIdx = RecvBuf_LBIdx + Recv_Disp[r];

      long *Match                  = new long [ SendEachRank_N[r] ];
      long *RecvPtr_LBIdx_IdxTable = new long [ SendEachRank_N[r] ];

//    be aware that there may be duplicate LBIdx in RecvPtr_LBIdx[]
//    --> RecvPtr_LBIdx_IdxTable[] is non-deterministic in that case
//    --> but it's not a real issue since the mapped PID is deterministic
      Mis_Heapsort( SendEachRank_N[r], RecvPtr_LBIdx, RecvPtr_LBIdx_IdxTable );

      Mis_Matching_int( (long)amr->NPatchComma[SonLv][1], amr->LB->IdxList_Real[SonLv], SendEachRank_N[r], RecvPtr_LBIdx, Match );

//    all target patches must be found
#     ifdef GAMER_DEBUG
      for (long t=0; t<SendEachRank_N[r]; t++)
         if ( Match[t] == -1L )
            Aux_Error( ERROR_INFO, "SonLv %d, TRank %d, LBIdx %ld found no matching patches !!\n",
                       SonLv, r, RecvPtr_LBIdx[t] );
#     endif


//    2.3.2. store the information for sending fine-grid B field
      SendEachRank_SibID[r] = new int [ SendEachRank_N[r] ];
      SendEachRank_PID  [r] = new int [ SendEachRank_N[r] ];

      for (long t=0; t<SendEachRank_N[r]; t++)
      {
         SendEachRank_SibID[r][t] = RecvPtr_SibID[t];

//       use index table since RecvPtr_LBIdx[] has been sorted
         const long sorted_idx = RecvPtr_LBIdx_IdxTable[t];
         SendEachRank_PID[r][sorted_idx]  = amr->LB->IdxList_Real_IdxTable[SonLv][ Match[t] ];
         SendEachRank_PID[r][sorted_idx] -= SendEachRank_PID[r][sorted_idx] % 8; // store the PID with LocalID=0
      }

      delete [] Match;
      delete [] RecvPtr_LBIdx_IdxTable;
   } // for (int r=0; r<MPI_NRank; r++)


// 2.4. free memory
   delete [] SendBuf_SibID;
   delete [] RecvBuf_SibID;
   delete [] SendBuf_LBIdx;
   delete [] RecvBuf_LBIdx;



// 3. transfer the fine-grid B field
// 3.1. set the MPI parameters
   real *SendBuf_FineB=NULL, *RecvBuf_FineB=NULL;

   for (int r=0; r<MPI_NRank; r++)
   {
       Send_NList[r] = SendEachRank_N[r]*(long)SQR( PS2 );
       Recv_NList[r] = RecvEachRank_N[r]*(long)SQR( PS2 );
   }

   Send_Disp[0] = 0L;
   Recv_Disp[0] = 0L;

   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp[r] = Send_Disp[r-1] + Send_NList[r-1];
      Recv_Disp[r] = Recv_Disp[r-1] + Recv_NList[r-1];
   }

   Send_NTotal = Send_Disp[ MPI_NRank-1 ] + Send_NList[ MPI_NRank-1 ];
   Recv_NTotal = Recv_Disp[ MPI_NRank-1 ] + Recv_NList[ MPI_NRank-1 ];

   SendBuf_FineB = new real [Send_NTotal];
   RecvBuf_FineB = new real [Recv_NTotal];


// 3.2. fill out the send buffer
   const int didx_in[6]     = { PS1, 0, SQR(PS1), 0, CUBE(PS1), 0 };    // x=PS1/0, y=PS1/0, z=PS1/0 faces
   const int TDir[3][2]     = { {1, 2}, {0, 2}, {0, 1} };               // transverse directions
   const int stride_in_n[3] = { PS1P1, 1, 1 };
   const int stride_in_m[3] = { PS1P1*PS1, PS1P1*PS1, PS1 };
   const int SonMagSg       = amr->MagSg[SonLv];

   for (int r=0; r<MPI_NRank; r++)
   {
      real *SendPtr0 = SendBuf_FineB + Send_Disp[r];

      for (long t=0; t<SendEachRank_N[r]; t++)
      {
         const int PID0 = SendEachRank_PID  [r][t];
         const int sib  = SendEachRank_SibID[r][t];
         const int dir  = sib/2;    // spatial direction

         real *SendPtr = SendPtr0 + t*SQR( PS2 );

//       check
#        ifdef GAMER_DEBUG
//       target patch must be a real patch
         if ( PID0 < 0  ||  PID0 >= amr->NPatchComma[SonLv][1] )
            Aux_Error( ERROR_INFO, "PID0 (%d) is not a real patch (SonLv %d, NRealPatch %d) !!\n",
                       PID0, SonLv, amr->NPatchComma[SonLv][1] );

//       target face must be a C-F interface
         const int sib_mirror = sib + 1 - 2*(sib&1);
         const int PID        = PID0 + TABLE_03( sib, 0 );
         const int SibPID     = amr->patch[0][SonLv][PID]->sibling[sib_mirror];
         if ( SibPID != -1 )
            Aux_Error( ERROR_INFO, "SibPID (%d) != -1 (SonLv %d, PID %d, sib_mirror %d) !!\n",
                       SibPID, SonLv, PID, sib_mirror );
#        endif

//       loop over the 4 sibling fine patches to collect the fine-grid B field on the C-F interfaces
         for (int g=0; g<4; g++)
         {
            const int LocalID    = TABLE_03( sib, g );
            const int PID        = PID0 + LocalID;
            const int didx_out_n = TABLE_02( LocalID, 'x'+TDir[dir][0], 0, PS1 );
            const int didx_out_m = TABLE_02( LocalID, 'x'+TDir[dir][1], 0, PS1 );

            for (int m=0; m<PS1; m++)  {  int idx_B_in  = m*stride_in_m[dir] + didx_in[sib];
                                          int idx_B_out = ( m + didx_out_m )*PS2 + didx_out_n;
            for (int n=0; n<PS1; n++)  {

               SendPtr[idx_B_out] = amr->patch[SonMagSg][SonLv][PID]->magnetic[dir][idx_B_in];

               idx_B_in  += stride_in_n[dir];
               idx_B_out ++;
            }}
         } // for (int g=0; g<4; g++)
      } // for (long t=0; t<SendEachRank_N[r]; t++)
   } // for (int r=0; r<MPI_NRank; r++)


// 3.3. invoke MPI
   MPI_Alltoallv_GAMER( SendBuf_FineB, Send_NList, Send_Disp, MPI_GAMER_REAL,
                        RecvBuf_FineB, Recv_NList, Recv_Disp, MPI_GAMER_REAL, MPI_COMM_WORLD );


// 3.4. set the data to be returned by this function
   CFB_BField = RecvBuf_FineB;

   for (int r=0; r<MPI_NRank; r++)  CFB_NSibEachRank[r] = RecvEachRank_N[r];


// 4. free memory
   delete [] SendBuf_FineB;

   for (int r=0; r<MPI_NRank; r++)
   {
      free( RecvEachRank_SibID    [r] );
      free( RecvEachRank_SibLBIdx [r] );
      delete [] SendEachRank_SibID[r];
      delete [] SendEachRank_PID  [r];
   }

} // FUNCTION : MHD_LB_Refine_GetCoarseFineInterfaceBField



#endif // #if ( MODEL == HYDRO  &&  defined MHD  &&  defined LOAD_BALANCE )
