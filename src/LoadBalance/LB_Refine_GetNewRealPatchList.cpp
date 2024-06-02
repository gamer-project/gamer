#include "GAMER.h"

#ifdef LOAD_BALANCE



void PrepareCData( const int FaLv, const int FaPID, real *const FaData,
                   const int FaSg_Flu, const int FaGhost_Flu, const int NSide_Flu,
                   const int FaSg_Pot, const int FaGhost_Pot, const int NSide_Pot,
                   const int FaSg_Mag, const int FaGhost_Mag,
                   const int BC_Face[], const int FluVarIdxList[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Refine_GetNewRealPatchList
// Description :  Get the lists of father patches at FaLv to allocate/deallocate real son patches at FaLv+1
//
// Note        :  1. This function is invoked by LB_Refine()
//                2. Coarse-grid data for creating new patches are also collected in NewCData_Away[]
//                   --> Data of all sibling-buffer patches at FaLv must be prepared in advance in order to
//                       prepare the coarse-grid data for spatial interpolation
//                3. Home/Away : target patches at home/not at home
//                4. Output New/DelCr1D_Away[] are sorted but NewCData_Away[] is unsorted
//                5. Use "call-by-reference" for the input parameters
//
// Parameter   :  FaLv                  : Target refinement level to be refined
//                NNew_Home             : Number of home patches at FaLv to allocate son patches
//                NewPID_Home           : Patch indices of home patches at FaLv to allocate son patches
//                NNew_Away             : Number of away patches at FaLv to allocate son patches
//                NewCr1D_Away          : Padded 1D corner of away patches at FaLv to allocate son patches
//                NewCr1D_Away_IdxTable : Index table of NewCr1D_Away[]
//                NewCData_Away         : Coarse-grid data of away patches at FaLv to allocate son patches
//                NDel_Home             : Number of home patches at FaLv to deallocate son patches
//                DelPID_Home           : Patch indices of home patches at FaLv to deallocate son patches
//                NDel_Away             : Number of away patches at FaLv to deallocate son patches
//                DelCr1D_Away          : Padded 1D corner of away patches at FaLv to deallocate son patches
//
//                PARTICLE-only parameters (call-by-reference)
//                RefineF2S_Send_NPatchTotal : Total number of patches for exchanging particles from fathers to sons
//                RefineF2S_Send_PIDList     : Patch indices for exchanging particles from fathers to sons
//
//                MHD-only parameters (call-by-reference)
//                CFB_SibLBIdx_Home : Load-balance indices of the siblings of home patches
//                CFB_SibLBIdx_Away : Load-balance indices of the siblings of away patches
//
//                Hybrid-scheme-only parameters (call-by-reference)
//                SwitchFinerLevelsToWaveScheme : Convert levels above level "FaLv"´ to wave scheme
//
// Return      :  NNew_Home, NewPID_Home, NNew_Away, NewCr1D_Away, NewCr1D_Away_IdxTable, NewCData_Away,
//                NDel_Home, DelPID_Home, NDel_Away, DelCr1D_Away, RefineF2S_Send_NPatchTotal, RefineF2S_Send_PIDList,
//                CFB_SibLBIdx_Home, CFB_SibLBIdx_Away
//-------------------------------------------------------------------------------------------------------
void LB_Refine_GetNewRealPatchList( const int FaLv, int &NNew_Home, int *&NewPID_Home, int &NNew_Away,
                                    ulong *&NewCr1D_Away, int *&NewCr1D_Away_IdxTable, real *&NewCData_Away,
                                    int &NDel_Home, int *&DelPID_Home, int &NDel_Away, ulong *&DelCr1D_Away,
                                    int &RefineF2S_Send_NPatchTotal, int *&RefineF2S_Send_PIDList,
                                    long (*&CFB_SibLBIdx_Home)[6], long (*&CFB_SibLBIdx_Away)[6],
                                    bool &SwitchFinerLevelsToWaveScheme )
{

// 1. construct the unsorted new/delete lists for real patches
// ==========================================================================================
   const int SonLv   = FaLv + 1;
   const int FaNReal = amr->NPatchComma[FaLv][1];
   const int MemUnit = 1 + FaNReal/MPI_NRank;         // set arbitrarily

   patch_t *TP=NULL;
   long   LBIdx, CP_Min_local, CP_Max_local;          // CP : CutPoint --> for resetting min/max of LB->CutPoint
   int    TRank;
   int    NewMemSize[MPI_NRank], NNew_Send[MPI_NRank];
   int    DelMemSize[MPI_NRank], NDel_Send[MPI_NRank];
   int   *NewPID_Send[MPI_NRank];
   ulong *NewCr1D_Send[MPI_NRank], *DelCr1D_Send[MPI_NRank];
#  ifdef MHD
   long (*CFB_SibLBIdx_Send[MPI_NRank])[6];
#  endif


// initialize variables
   for (int r=0; r<MPI_NRank; r++)
   {
      NewMemSize       [r] = MemUnit;
      DelMemSize       [r] = MemUnit;
      NewCr1D_Send     [r] = ( ulong*     )malloc( NewMemSize[r]*sizeof(ulong)   );
      DelCr1D_Send     [r] = ( ulong*     )malloc( DelMemSize[r]*sizeof(ulong)   );
      NewPID_Send      [r] = ( int*       )malloc( NewMemSize[r]*sizeof(int  )   );
#     ifdef MHD
      CFB_SibLBIdx_Send[r] = ( long(*)[6] )malloc( NewMemSize[r]*sizeof(long )*6 );
#     endif
      NNew_Send        [r] = 0;
      NDel_Send        [r] = 0;
   }

   NewPID_Home       = ( int*       )malloc( NewMemSize[MPI_Rank]*sizeof(int )   );
   DelPID_Home       = ( int*       )malloc( DelMemSize[MPI_Rank]*sizeof(int )   );
#  ifdef MHD
   CFB_SibLBIdx_Home = ( long(*)[6] )malloc( NewMemSize[MPI_Rank]*sizeof(long)*6 );
#  endif
   NNew_Home         = 0;
   NDel_Home         = 0;


// loop over all real patches at FaLv
   for (int FaPID=0; FaPID<FaNReal; FaPID++)
   {
      TP = amr->patch[0][FaLv][FaPID];

//    1.1 allocate list
      if ( TP->flag  &&  TP->son == -1 )
      {
//       set the target rank
//###NOTE: faster version can only be applied to the Hilbert space-filling curve
#        if ( LOAD_BALANCE == HILBERT )
         LBIdx = 8*TP->LB_Idx;   // faster
#        else
         LBIdx = LB_Corner2Index( SonLv, TP->corner, CHECK_ON );
#        endif
         TRank = LB_Index2Rank( SonLv, LBIdx, CHECK_OFF );

//       LB_Idx of the newly-created sons can lie outside LB->CutPoint --> need to reset LB->CutPoint
         if ( TRank == -1 )
         {
            if ( LBIdx < amr->LB->CutPoint[SonLv][0] )
            {
               TRank = 0;
               amr->LB->CutPoint[SonLv][0] = LBIdx - LBIdx%8;
            }

            else if ( LBIdx >= amr->LB->CutPoint[SonLv][MPI_NRank] )
            {
               TRank = MPI_NRank - 1;
               amr->LB->CutPoint[SonLv][MPI_NRank] = LBIdx - LBIdx%8 + 8;
            }

#           ifdef GAMER_DEBUG
            else
               Aux_Error( ERROR_INFO, "SonLv %d, incorrect LBIdx %ld, LBIdx_Min %ld, LBIdx_Max %ld !!\n",
                          SonLv, LBIdx, amr->LB->CutPoint[SonLv][0], amr->LB->CutPoint[SonLv][MPI_NRank] );
#           endif
         }


//       get the sibling information required for the B field interpolation in MHD
#        ifdef MHD
         long CFB_SibLBIdx[6] = { -1, -1, -1, -1, -1, -1 };

         for (int s=0; s<6; s++)
         {
            const int SibPID    = amr->patch[0][FaLv][FaPID]->sibling[s];
            const int SibSonPID = ( SibPID >= 0 ) ? amr->patch[0][FaLv][SibPID]->son : -1;

            if ( SibSonPID != -1 )
            {
#              if ( LOAD_BALANCE == HILBERT )
               CFB_SibLBIdx[s] = 8*amr->patch[0][FaLv][SibPID]->LB_Idx; // faster
#              else
               CFB_SibLBIdx[s] = LB_Corner2Index( SonLv, amr->patch[0][FaLv][SibPID]->corner, CHECK_ON );
#              endif
            }
         }
#        endif // #ifdef MHD

#        if ( ELBDM_SCHEME == ELBDM_HYBRID )
         if ( !amr->use_wave_flag[SonLv]  &&  TP->switch_to_wave_flag ) {
            SwitchFinerLevelsToWaveScheme = true;
         }
#        endif


//       record the new lists
         if ( TRank == MPI_Rank ) // target son patches are home
         {
//          allocate enough memory
            if ( NNew_Home >= NewMemSize[TRank] )
            {
               NewMemSize[TRank] += MemUnit;
               NewPID_Home        = ( int*       )realloc( NewPID_Home,       NewMemSize[TRank]*sizeof(int )   );
#              ifdef MHD
               CFB_SibLBIdx_Home  = ( long(*)[6] )realloc( CFB_SibLBIdx_Home, NewMemSize[TRank]*sizeof(long)*6 );
#              endif
            }

            NewPID_Home      [NNew_Home]    = FaPID;
#           ifdef MHD
            for (int s=0; s<6; s++)
            CFB_SibLBIdx_Home[NNew_Home][s] = CFB_SibLBIdx[s];
#           endif

            NNew_Home ++;
         } // if ( TRank == MPI_Rank )

         else // TRank != MPI_Rank (target son patches are not home)
         {
//          allocate enough memory
            if ( NNew_Send[TRank] >= NewMemSize[TRank] )
            {
               NewMemSize       [TRank] += MemUnit;
               NewCr1D_Send     [TRank]  = ( ulong*     )realloc( NewCr1D_Send     [TRank], NewMemSize[TRank]*sizeof(ulong)   );
               NewPID_Send      [TRank]  = ( int*       )realloc( NewPID_Send      [TRank], NewMemSize[TRank]*sizeof(int  )   );
#              ifdef MHD
               CFB_SibLBIdx_Send[TRank]  = ( long(*)[6] )realloc( CFB_SibLBIdx_Send[TRank], NewMemSize[TRank]*sizeof(long )*6 );
#              endif
            }

            NewCr1D_Send     [TRank][ NNew_Send[TRank] ]    = TP->PaddedCr1D;
            NewPID_Send      [TRank][ NNew_Send[TRank] ]    = FaPID;
#           ifdef MHD
            for (int s=0; s<6; s++)
            CFB_SibLBIdx_Send[TRank][ NNew_Send[TRank] ][s] = CFB_SibLBIdx[s];
#           endif

            NNew_Send[TRank] ++;

#           ifdef PARTICLE
#           ifdef DEBUG_PARTICLE
            if ( RefineF2S_Send_NPatchTotal >= FaNReal )
               Aux_Error( ERROR_INFO, "target index (%d) >= FaNReal (%d) !!\n", RefineF2S_Send_NPatchTotal, FaNReal );
#           endif

            RefineF2S_Send_PIDList[ RefineF2S_Send_NPatchTotal ++ ] = FaPID;
#           endif // #ifdef PARTICLE
         } // if ( TRank == MPI_Rank ) ... else ...
      } // if ( TP->flag  &&  TP->son == -1 )


//    1.2 deallocate list
      else if ( !TP->flag  &&  TP->son != -1 )
      {
//       set the target rank
//###NOTE: faster version can only be applied to the Hilbert space-filling curve
#        if ( LOAD_BALANCE == HILBERT )
         LBIdx = 8*TP->LB_Idx;   // faster
#        else
         LBIdx = LB_Corner2Index( SonLv, TP->corner, CHECK_ON );
#        endif
         TRank = LB_Index2Rank( SonLv, LBIdx, CHECK_ON );


//       record the delete list
         if ( TRank == MPI_Rank ) // target son patches are home
         {
//          allocate enough memory
            if ( NDel_Home >= DelMemSize[TRank] )
            {
               DelMemSize[TRank] += MemUnit;
               DelPID_Home        = (int*)realloc( DelPID_Home, DelMemSize[TRank]*sizeof(int) );
            }

            DelPID_Home[ NDel_Home ++ ] = FaPID;

//          check : SonPID should be home
#           ifdef GAMER_DEBUG
            if ( TP->son < -1 )
               Aux_Error( ERROR_INFO, "FaLv %d, FaPID %d, SonPID = %d is not home !!\n", FaLv, FaPID, TP->son );
#           endif
         }

         else // TRank != MPI_Rank (target son patches are not home)
         {
//          allocate enough memory
            if ( NDel_Send[TRank] >= DelMemSize[TRank] )
            {
               DelMemSize  [TRank] += MemUnit;
               DelCr1D_Send[TRank]  = (ulong*)realloc( DelCr1D_Send[TRank], DelMemSize[TRank]*sizeof(ulong) );
            }

            DelCr1D_Send[TRank][ NDel_Send[TRank] ++ ] = TP->PaddedCr1D;

//          check : SonPID should not be home
#           ifdef GAMER_DEBUG
            if ( TP->son >= 0 )
               Aux_Error( ERROR_INFO, "FaLv %d, FaPID %d, SonPID = %d is home !!\n", FaLv, FaPID, TP->son );
#           endif
         }

      } // else if ( !TP->flag  &&  TP->son != -1 )
   } // for (int FaPID=0; FaPID<FaNReal; FaPID++)


// 1.3 get the new LB_CutPoint at SonLv
   CP_Min_local = amr->LB->CutPoint[SonLv][        0];
   CP_Max_local = amr->LB->CutPoint[SonLv][MPI_NRank];

   MPI_Allreduce( &CP_Min_local, &amr->LB->CutPoint[SonLv][        0], 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &CP_Max_local, &amr->LB->CutPoint[SonLv][MPI_NRank], 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD );



// 2. broadcast the unsorted new/delete lists to all other ranks
// ============================================================================================================
   const int FaSg_Flu = amr->FluSg[FaLv];
   int NSide_Flu, FaGhost_Flu, FaSize_Flu;
   int PSize=0;

   Int_Table( OPT__REF_FLU_INT_SCHEME, NSide_Flu, FaGhost_Flu );

   FaSize_Flu = PS1 + 2*FaGhost_Flu;
   PSize     += NCOMP_TOTAL*CUBE( FaSize_Flu );

#  ifdef GRAVITY
   const int FaSg_Pot = amr->PotSg[FaLv];
   int NSide_Pot, FaGhost_Pot, FaSize_Pot;

   Int_Table( OPT__REF_POT_INT_SCHEME, NSide_Pot, FaGhost_Pot );

   FaSize_Pot = PS1 + 2*FaGhost_Pot;
   PSize     += CUBE( FaSize_Pot );
#  else
   const int FaSg_Pot=NULL_INT, NSide_Pot=NULL_INT, FaGhost_Pot=NULL_INT, FaSize_Pot=NULL_INT;
#  endif

#  ifdef MHD
   const int FaSg_Mag = amr->MagSg[FaLv];
   int NSide_Mag_Useless, FaGhost_Mag, FaSize_Mag_T, FaSize_Mag_N;

   Int_Table( OPT__REF_MAG_INT_SCHEME, NSide_Mag_Useless, FaGhost_Mag );

   FaSize_Mag_T = PS1 + 2*FaGhost_Mag;    // coarse-grid size along the transverse (_T) / normal (_N) direction
   FaSize_Mag_N = PS1P1;
   PSize       += NCOMP_MAG*FaSize_Mag_N*SQR( FaSize_Mag_T );
#  else
   const int FaSg_Mag=NULL_INT, FaGhost_Mag=NULL_INT, FaSize_Mag_T=NULL_INT, FaSize_Mag_N=NULL_INT;
#  endif


   int    New_Send_Disp[MPI_NRank], New_Recv_Disp[MPI_NRank], NNew_Recv[MPI_NRank], NNew_Send_Total, NNew_Recv_Total;
   int    Del_Send_Disp[MPI_NRank], Del_Recv_Disp[MPI_NRank], NDel_Recv[MPI_NRank], NDel_Send_Total, NDel_Recv_Total;
   int    Counter;
   long   NNew_Send_CData[MPI_NRank], NNew_Recv_CData[MPI_NRank];
   long   New_Send_Disp_CData[MPI_NRank], New_Recv_Disp_CData[MPI_NRank];
   ulong *New_SendBuf_Cr1D=NULL, *New_RecvBuf_Cr1D=NULL, *Del_SendBuf_Cr1D=NULL, *Del_RecvBuf_Cr1D=NULL;
   real  *New_SendBuf_CData=NULL, *New_RecvBuf_CData=NULL;

#  ifdef MHD
   int    CFB_Send_Disp_SibLBIdx[MPI_NRank], CFB_Recv_Disp_SibLBIdx[MPI_NRank];
   int    CFB_Send_NList_SibLBIdx[MPI_NRank], CFB_Recv_NList_SibLBIdx[MPI_NRank];
   long (*CFB_SendBuf_SibLBIdx)[6]=NULL, (*CFB_RecvBuf_SibLBIdx)[6]=NULL;
#  endif

// 2.1 broadcast the number of elements sent to different ranks
   MPI_Alltoall( NNew_Send, 1, MPI_INT, NNew_Recv, 1, MPI_INT, MPI_COMM_WORLD );
   MPI_Alltoall( NDel_Send, 1, MPI_INT, NDel_Recv, 1, MPI_INT, MPI_COMM_WORLD );


// 2.2 allocate the MPI buffers
   New_Send_Disp[0] = 0;
   New_Recv_Disp[0] = 0;
   Del_Send_Disp[0] = 0;
   Del_Recv_Disp[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      New_Send_Disp[r] = New_Send_Disp[r-1] + NNew_Send[r-1];
      New_Recv_Disp[r] = New_Recv_Disp[r-1] + NNew_Recv[r-1];
      Del_Send_Disp[r] = Del_Send_Disp[r-1] + NDel_Send[r-1];
      Del_Recv_Disp[r] = Del_Recv_Disp[r-1] + NDel_Recv[r-1];
   }

   NNew_Send_Total = New_Send_Disp[MPI_NRank-1] + NNew_Send[MPI_NRank-1];
   NNew_Recv_Total = New_Recv_Disp[MPI_NRank-1] + NNew_Recv[MPI_NRank-1];
   NDel_Send_Total = Del_Send_Disp[MPI_NRank-1] + NDel_Send[MPI_NRank-1];
   NDel_Recv_Total = Del_Recv_Disp[MPI_NRank-1] + NDel_Recv[MPI_NRank-1];

   for (int r=0; r<MPI_NRank; r++)
   {
      NNew_Send_CData        [r] = (long)PSize*(long)NNew_Send    [r];
      NNew_Recv_CData        [r] = (long)PSize*(long)NNew_Recv    [r];
      New_Send_Disp_CData    [r] = (long)PSize*(long)New_Send_Disp[r];
      New_Recv_Disp_CData    [r] = (long)PSize*(long)New_Recv_Disp[r];
#     ifdef MHD
      CFB_Send_NList_SibLBIdx[r] =     6*NNew_Send    [r];
      CFB_Recv_NList_SibLBIdx[r] =     6*NNew_Recv    [r];
      CFB_Send_Disp_SibLBIdx [r] =     6*New_Send_Disp[r];
      CFB_Recv_Disp_SibLBIdx [r] =     6*New_Recv_Disp[r];
#     endif
   }

// variables to be returned by this function
   NNew_Away             = NNew_Recv_Total;
   NDel_Away             = NDel_Recv_Total;
   NewCr1D_Away          = new ulong [NNew_Recv_Total      ];
   NewCr1D_Away_IdxTable = new int   [NNew_Recv_Total      ];
   NewCData_Away         = new real  [NNew_Recv_Total*PSize];
   DelCr1D_Away          = new ulong [NDel_Recv_Total      ];
   New_SendBuf_Cr1D      = new ulong [NNew_Send_Total      ];
   New_RecvBuf_Cr1D      = NewCr1D_Away;
   New_SendBuf_CData     = new real  [NNew_Send_Total*PSize];
   New_RecvBuf_CData     = NewCData_Away;
   Del_SendBuf_Cr1D      = new ulong [NDel_Send_Total      ];
   Del_RecvBuf_Cr1D      = DelCr1D_Away;
#  ifdef MHD
   CFB_SendBuf_SibLBIdx  = new long  [NNew_Send_Total      ][6];
   CFB_RecvBuf_SibLBIdx  = new long  [NNew_Recv_Total      ][6];
#  endif


// 2.3 prepare the MPI send buffers
// 2.3.1&2 new Cr1D/CData (and SibLBIdx for MHD)

// determine the priority of different boundary faces (z>y>x) to set the corner cells properly for the non-periodic B.C.
   int BC_Face[26], BC_Face_tmp[3], FluVarIdxList[NCOMP_TOTAL];

   for (int s=0; s<26; s++)
   {
      BC_Face_tmp[0] = TABLE_01( s, 'x', 0, -1, 1 );
      BC_Face_tmp[1] = TABLE_01( s, 'y', 2, -1, 3 );
      BC_Face_tmp[2] = TABLE_01( s, 'z', 4, -1, 5 );

//    z > y > x
      if      ( BC_Face_tmp[2] != -1  &&  OPT__BC_FLU[BC_Face_tmp[2]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[2];
      else if ( BC_Face_tmp[1] != -1  &&  OPT__BC_FLU[BC_Face_tmp[1]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[1];
      else if ( BC_Face_tmp[0] != -1  &&  OPT__BC_FLU[BC_Face_tmp[0]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[0];
      else                                                                                   BC_Face[s] = NULL_INT;
   }

   for (int v=0; v<NCOMP_TOTAL; v++)   FluVarIdxList[v] = v;

// prepare the coarse-grid data
   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<NNew_Send[r]; t++)
   {
      New_SendBuf_Cr1D    [Counter]    = NewCr1D_Send    [r][t];
#     ifdef MHD
      for (int s=0; s<6; s++)
      CFB_SendBuf_SibLBIdx[Counter][s] = CFB_SibLBIdx_Send[r][t][s];
#     endif

      PrepareCData( FaLv, NewPID_Send[r][t], New_SendBuf_CData+Counter*PSize,
                    FaSg_Flu, FaGhost_Flu, NSide_Flu, FaSg_Pot, FaGhost_Pot, NSide_Pot, FaSg_Mag, FaGhost_Mag,
                    BC_Face, FluVarIdxList );

      Counter ++;
   }


// 2.3.3 delete Cr1D
   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<NDel_Send[r]; t++)
      Del_SendBuf_Cr1D[ Counter ++ ] = DelCr1D_Send[r][t];


// 2.4 broadcast the send data
// 2.4.1 new Cr1D
   MPI_Alltoallv( New_SendBuf_Cr1D, NNew_Send, New_Send_Disp, MPI_UNSIGNED_LONG,
                  New_RecvBuf_Cr1D, NNew_Recv, New_Recv_Disp, MPI_UNSIGNED_LONG, MPI_COMM_WORLD );

// 2.4.2 SibLBIdx for MHD
#  ifdef MHD
   MPI_Alltoallv( CFB_SendBuf_SibLBIdx, CFB_Send_NList_SibLBIdx, CFB_Send_Disp_SibLBIdx, MPI_LONG,
                  CFB_RecvBuf_SibLBIdx, CFB_Recv_NList_SibLBIdx, CFB_Recv_Disp_SibLBIdx, MPI_LONG, MPI_COMM_WORLD );
#  endif

// 2.4.3 new CData
   MPI_Alltoallv_GAMER( New_SendBuf_CData, NNew_Send_CData, New_Send_Disp_CData, MPI_GAMER_REAL,
                        New_RecvBuf_CData, NNew_Recv_CData, New_Recv_Disp_CData, MPI_GAMER_REAL, MPI_COMM_WORLD );

// 2.4.4 delete Cr1D
   MPI_Alltoallv( Del_SendBuf_Cr1D, NDel_Send, Del_Send_Disp, MPI_UNSIGNED_LONG,
                  Del_RecvBuf_Cr1D, NDel_Recv, Del_Recv_Disp, MPI_UNSIGNED_LONG, MPI_COMM_WORLD );



// 3. sort *Cr1D_Away[] and CFB_SibLBIdx_Away[]
// ============================================================================================================
   Mis_Heapsort           ( NNew_Away, NewCr1D_Away, NewCr1D_Away_IdxTable );
   Mis_Heapsort<int,ulong>( NDel_Away, DelCr1D_Away, NULL                  );

// must ensure that CFB_SibLBIdx_Away[] and Cr1D_Away[] are sorted consistently
// --> so that the coarse-fine interface B field data received in MHD_LB_Refine_GetCoarseFineInterfaceBField()
//     are in the correct order
#  ifdef MHD
   const long (*CFB_SibLBIdx_Away_Unsorted)[6] = CFB_RecvBuf_SibLBIdx;
   CFB_SibLBIdx_Away = new long [NNew_Away][6];

   for (int idx_sorted=0; idx_sorted<NNew_Away; idx_sorted++)
   {
      const int idx_unsorted = NewCr1D_Away_IdxTable[idx_sorted];

      for (int s=0; s<6; s++)
         CFB_SibLBIdx_Away[idx_sorted][s] = CFB_SibLBIdx_Away_Unsorted[idx_unsorted][s];
   }
#  endif



// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( NewCr1D_Send    [r] );
      free( DelCr1D_Send    [r] );
      free( NewPID_Send     [r] );
#     ifdef MHD
      free( CFB_SibLBIdx_Send[r] );
#     endif
   }
   delete [] New_SendBuf_Cr1D;
   delete [] New_SendBuf_CData;
   delete [] Del_SendBuf_Cr1D;
#  ifdef MHD
   delete [] CFB_SendBuf_SibLBIdx;
   delete [] CFB_RecvBuf_SibLBIdx;
#  endif

} // FUNCTION : LB_Refine_GetNewRealPatchList



//-------------------------------------------------------------------------------------------------------
// Function    :  PrepareCData
// Description :  Prepare coarse-grid data for spatial interpolation
//
// Note        :  1. Data of all sibling-buffer patches at FaLv must be prepared in advance
//                2. This function is also used by LB_Refine_AllocateNewPatch()
//
// Parameter   :  FaLv          : Coarse-grid refinement level
//                FaPID         : Father patch index to prepare the coarse-grid data
//                FaData        : Array to store the coarse-grid data
//                FaSg_Flu      : Sandglass of the fluid data
//                FaGhost_Flu   : Ghost size of the fluid data
//                NSide_Flu     : Number of sibling directions to prepare the ghost-zone data (6/26) for the fluid data
//                FaSg_Pot      : Sandglass of the potential data
//                FaGhost_Pot   : Ghost size of the potential data
//                NSide_Pot     : Number of sibling directions to prepare the ghost-zone data (6/26) for the potential data
//                FaSg_Mag      : Sandglass of the magnetic field
//                FaGhost_Mag   : Ghost size of the magnetic field (only for the transverse direction)
//                BC_Face       : Corresponding boundary faces (0~5) along 26 sibling directions -> for non-periodic B.C. only
//                FluVarIdxList : List of target fluid variable indices                          -> for non-periodic B.C. only
//
// Return      :  FaData
//-------------------------------------------------------------------------------------------------------
void PrepareCData( const int FaLv, const int FaPID, real *const FaData,
                   const int FaSg_Flu, const int FaGhost_Flu, const int NSide_Flu,
                   const int FaSg_Pot, const int FaGhost_Pot, const int NSide_Pot,
                   const int FaSg_Mag, const int FaGhost_Mag,
                   const int BC_Face[], const int FluVarIdxList[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( Flu_ParaBuf < FaGhost_Flu )
      Aux_Error( ERROR_INFO, "Flu_ParaBuf (%d) < FaGhost_Flu (%d) --> refinement will fail !!\n",
                 Flu_ParaBuf, FaGhost_Flu );
#  ifdef GRAVITY
   if ( Pot_ParaBuf < FaGhost_Pot )
      Aux_Error( ERROR_INFO, "Pot_ParaBuf (%d) < FaGhost_Pot (%d) --> refinement will fail !!\n",
                 Pot_ParaBuf, FaGhost_Pot );
#  endif
#  ifdef MHD
// Flu_ParaBuf is used for transferring B field as well
   if ( Flu_ParaBuf < FaGhost_Mag )
      Aux_Error( ERROR_INFO, "Flu_ParaBuf (%d) < FaGhost_Mag (%d) --> refinement will fail !!\n",
                 Flu_ParaBuf, FaGhost_Mag );
#  endif
#  endif // #ifdef GAMER_DEBUG


// 1. fill up the central region of FaData
   real *FaData_Next = FaData;

   const int FaSize_Flu    = PS1 + 2*FaGhost_Flu;
   real *const FaData_Flu  = FaData_Next;
   FaData_Next            += NCOMP_TOTAL*CUBE( FaSize_Flu );

#  ifdef GRAVITY
   const int FaSize_Pot    = PS1 + 2*FaGhost_Pot;
   real *const FaData_Pot  = FaData_Next;
   FaData_Next            += CUBE( FaSize_Pot );
#  endif

#  ifdef MHD
   const int FaSize_Mag_T  = PS1 + 2*FaGhost_Mag;  // coarse-grid size along the transverse (_T) / normal (_N) direction
   const int FaSize_Mag_N  = PS1P1;
   real *const FaData_MagX = FaData_Next + MAGX*FaSize_Mag_N*SQR( FaSize_Mag_T );
   real *const FaData_MagY = FaData_Next + MAGY*FaSize_Mag_N*SQR( FaSize_Mag_T );
   real *const FaData_MagZ = FaData_Next + MAGZ*FaSize_Mag_N*SQR( FaSize_Mag_T );
   FaData_Next            += NCOMP_MAG*FaSize_Mag_N*SQR( FaSize_Mag_T );
#  endif

   int idx_out, i_out, j_out, k_out;


// 1.1 fluid data
   for (int v=0; v<NCOMP_TOTAL; v++)  {
   for (int k=0; k<PS1; k++)  {  k_out = k + FaGhost_Flu;
   for (int j=0; j<PS1; j++)  {  j_out = j + FaGhost_Flu;
   for (int i=0; i<PS1; i++)  {  i_out = i + FaGhost_Flu;

      idx_out = ((v*FaSize_Flu + k_out)*FaSize_Flu + j_out)*FaSize_Flu + i_out;

      FaData_Flu[idx_out] = amr->patch[FaSg_Flu][FaLv][FaPID]->fluid[v][k][j][i];

   }}}}


// 1.2 potential data
#  ifdef GRAVITY
   for (int k=0; k<PS1; k++)  {  k_out = k + FaGhost_Pot;
   for (int j=0; j<PS1; j++)  {  j_out = j + FaGhost_Pot;
   for (int i=0; i<PS1; i++)  {  i_out = i + FaGhost_Pot;

      idx_out = (k_out*FaSize_Pot + j_out)*FaSize_Pot + i_out;

      FaData_Pot[idx_out] = amr->patch[FaSg_Pot][FaLv][FaPID]->pot[k][j][i];

   }}}
#  endif


// 1.3 magnetic field
#  ifdef MHD
   int idx_B_in, idx_B_out;

// Bx
   idx_B_in = 0;
   for (int k=FaGhost_Mag; k<FaGhost_Mag+PS1; k++)  {
   for (int j=FaGhost_Mag; j<FaGhost_Mag+PS1; j++)  {  idx_B_out = IDX321(           0, j, k, FaSize_Mag_N, FaSize_Mag_T );
   for (int i=0;           i<FaSize_Mag_N;    i++)  {
      FaData_MagX[ idx_B_out ++ ] = amr->patch[FaSg_Mag][FaLv][FaPID]->magnetic[MAGX][ idx_B_in ++ ];
   }}}

// By
   idx_B_in = 0;
   for (int k=FaGhost_Mag; k<FaGhost_Mag+PS1; k++)  {
   for (int j=0;           j<FaSize_Mag_N;    j++)  {  idx_B_out = IDX321( FaGhost_Mag, j, k, FaSize_Mag_T, FaSize_Mag_N );
   for (int i=FaGhost_Mag; i<FaGhost_Mag+PS1; i++)  {
      FaData_MagY[ idx_B_out ++ ] = amr->patch[FaSg_Mag][FaLv][FaPID]->magnetic[MAGY][ idx_B_in ++ ];
   }}}

// Bz
   idx_B_in = 0;
   for (int k=0;           k<FaSize_Mag_N;    k++)  {
   for (int j=FaGhost_Mag; j<FaGhost_Mag+PS1; j++)  {  idx_B_out = IDX321( FaGhost_Mag, j, k, FaSize_Mag_T, FaSize_Mag_T );
   for (int i=FaGhost_Mag; i<FaGhost_Mag+PS1; i++)  {
      FaData_MagZ[ idx_B_out ++ ] = amr->patch[FaSg_Mag][FaLv][FaPID]->magnetic[MAGZ][ idx_B_in ++ ];
   }}}
#  endif // #ifdef MHD



// 2. fill up the ghost zones of FaData (no interpolation is required)
   const int   NDer       = 0;
   const long *DerVarList = NULL;

   int    loop[3], offset_out[3], offset_in[3], i_in, j_in, k_in;
   int    BC_Sibling, BC_Idx_Start[3], BC_Idx_End[3];
   double xyz_flu[3];

// calculate the corner coordinates of the coarse-grid data for the user-specified B.C.
   for (int d=0; d<3; d++)    xyz_flu[d] = amr->patch[0][FaLv][FaPID]->EdgeL[d] + (0.5-FaGhost_Flu)*amr->dh[FaLv];


// 2.1 fluid data
   for (int sib=0; sib<NSide_Flu; sib++)
   {
      const int SibPID = amr->patch[0][FaLv][FaPID]->sibling[sib];

      for (int d=0; d<3; d++)
      {
         loop      [d] = TABLE_01( sib, 'x'+d, FaGhost_Flu, PS1, FaGhost_Flu );
         offset_out[d] = TABLE_01( sib, 'x'+d, 0, FaGhost_Flu, FaGhost_Flu+PS1 );
      }

//    2.1.1 if the target sibling patch exists --> just copy data from it directly
      if ( SibPID >= 0 )
      {
         for (int d=0; d<3; d++)    offset_in[d] = TABLE_01( sib, 'x'+d, PS1-FaGhost_Flu, 0, 0 );

         for (int v=0; v<NCOMP_TOTAL; v++)  {
         for (int k=0; k<loop[2]; k++)  {  k_out = k + offset_out[2];  k_in = k + offset_in[2];
         for (int j=0; j<loop[1]; j++)  {  j_out = j + offset_out[1];  j_in = j + offset_in[1];
         for (int i=0; i<loop[0]; i++)  {  i_out = i + offset_out[0];  i_in = i + offset_in[0];

            idx_out = ((v*FaSize_Flu + k_out)*FaSize_Flu + j_out)*FaSize_Flu + i_out;

            FaData_Flu[idx_out] = amr->patch[FaSg_Flu][FaLv][SibPID]->fluid[v][k_in][j_in][i_in];

         }}}}
      }


//    2.1.2 if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
      else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
      {
         for (int d=0; d<3; d++)
         {
            BC_Idx_Start[d] = offset_out[d];
            BC_Idx_End  [d] = loop[d] + BC_Idx_Start[d] - 1;
         }

         BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;

#        ifdef GAMER_DEBUG
         if ( BC_Face[BC_Sibling] < 0  ||  BC_Face[BC_Sibling] > 5 )
            Aux_Error( ERROR_INFO, "incorrect BC_Face[%d] = %d !!\n", BC_Sibling, BC_Face[BC_Sibling] );

         if ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] == BC_FLU_PERIODIC )
            Aux_Error( ERROR_INFO, "OPT__BC_FLU == BC_FLU_PERIODIC (BC_Sibling %d, BC_Face %d, SibPID %d, FaPID %d, sib %d, FaLv %d) !!\n",
                       BC_Sibling, BC_Face[BC_Sibling], SibPID, FaPID, sib, FaLv );
#        endif

         switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
         {
#           if ( MODEL == HYDRO )
            case BC_FLU_OUTFLOW:
               Hydro_BoundaryCondition_Outflow   ( FaData_Flu, BC_Face[BC_Sibling], NCOMP_TOTAL, FaGhost_Flu,
                                                   FaSize_Flu, FaSize_Flu, FaSize_Flu, BC_Idx_Start, BC_Idx_End );
            break;

            case BC_FLU_REFLECTING:
               Hydro_BoundaryCondition_Reflecting( FaData_Flu, BC_Face[BC_Sibling], NCOMP_TOTAL, FaGhost_Flu,
                                                   FaSize_Flu, FaSize_Flu, FaSize_Flu, BC_Idx_Start, BC_Idx_End,
                                                   FluVarIdxList, NDer, DerVarList );
            break;

            case BC_FLU_DIODE:
               Hydro_BoundaryCondition_Diode(      FaData_Flu, BC_Face[BC_Sibling], NCOMP_TOTAL, FaGhost_Flu,
                                                   FaSize_Flu, FaSize_Flu, FaSize_Flu, BC_Idx_Start, BC_Idx_End,
                                                   FluVarIdxList, NDer, DerVarList );
            break;
#           endif

            case BC_FLU_USER:
               Flu_BoundaryCondition_User        ( FaData_Flu,                      NCOMP_TOTAL, FaGhost_Flu,
                                                   FaSize_Flu, FaSize_Flu, FaSize_Flu, BC_Idx_Start, BC_Idx_End,
                                                   FluVarIdxList, Time[FaLv], amr->dh[FaLv], xyz_flu, _TOTAL, FaLv );
            break;

            default:
               Aux_Error( ERROR_INFO, "unsupported fluid B.C. (%d) !!\n", OPT__BC_FLU[ BC_Face[BC_Sibling] ] );

         } // switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
      } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )


//    2.1.3 it will violate the proper-nesting condition if the flagged patch is NOT surrounded by siblings
      else
         Aux_Error( ERROR_INFO, "SibPID == %d (FaLv %d, FaPID %d, sib %d) !!\n", SibPID, FaLv, FaPID, sib );

   } // for (int sib=0; sib<NSide_Flu; sib++)


// 2.2 potential data
#  ifdef GRAVITY
   for (int sib=0; sib<NSide_Pot; sib++)
   {
      const int SibPID = amr->patch[0][FaLv][FaPID]->sibling[sib];

      for (int d=0; d<3; d++)
      {
         loop      [d] = TABLE_01( sib, 'x'+d, FaGhost_Pot, PS1, FaGhost_Pot );
         offset_out[d] = TABLE_01( sib, 'x'+d, 0, FaGhost_Pot, FaGhost_Pot+PS1 );
      }

//    2.2.1 if the target sibling patch exists --> just copy data from it directly
      if ( SibPID >= 0 )
      {
         for (int d=0; d<3; d++)    offset_in[d] = TABLE_01( sib, 'x'+d, PS1-FaGhost_Pot, 0, 0 );

         for (int k=0; k<loop[2]; k++)  {  k_out = k + offset_out[2];  k_in = k + offset_in[2];
         for (int j=0; j<loop[1]; j++)  {  j_out = j + offset_out[1];  j_in = j + offset_in[1];
         for (int i=0; i<loop[0]; i++)  {  i_out = i + offset_out[0];  i_in = i + offset_in[0];

            idx_out = (k_out*FaSize_Pot + j_out)*FaSize_Pot + i_out;

            FaData_Pot[idx_out] = amr->patch[FaSg_Pot][FaLv][SibPID]->pot[k_in][j_in][i_in];

         }}}
      }


//    2.2.2 if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
      else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
      {
         for (int d=0; d<3; d++)
         {
            BC_Idx_Start[d] = offset_out[d];
            BC_Idx_End  [d] = loop[d] + BC_Idx_Start[d] - 1;
         }

         BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;

#        ifdef GAMER_DEBUG
         if ( BC_Face[BC_Sibling] < 0  ||  BC_Face[BC_Sibling] > 5 )
            Aux_Error( ERROR_INFO, "incorrect BC_Face[%d] = %d !!\n", BC_Sibling, BC_Face[BC_Sibling] );

         if ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] == BC_FLU_PERIODIC )
            Aux_Error( ERROR_INFO, "OPT__BC_FLU == BC_FLU_PERIODIC (BC_Sibling %d, BC_Face %d, SibPID %d, FaPID %d, sib %d, FaLv %d) !!\n",
                       BC_Sibling, BC_Face[BC_Sibling], SibPID, FaPID, sib, FaLv );
#        endif

//       extrapolate potential
         Poi_BoundaryCondition_Extrapolation( FaData_Pot, BC_Face[BC_Sibling], 1, FaGhost_Pot,
                                              FaSize_Pot, FaSize_Pot, FaSize_Pot, BC_Idx_Start, BC_Idx_End );
      }


//    2.2.3 it will violate the proper-nesting condition if the flagged patch is NOT surrounded by siblings
      else
         Aux_Error( ERROR_INFO, "SibPID == %d (FaLv %d, FaPID %d, sib %d) !!\n", SibPID, FaLv, FaPID, sib );

   } // for (int sib=0; sib<NSide_Pot; sib++)
#  endif // #ifdef GRAVITY


// 2.3 magnetic field
#  ifdef MHD
// interpolation on B field only requires ghost zones along the two transverse directions
// --> skip sib>=6 since ghost zones along the diagonal directions are not required
   for (int sib=0; sib<6; sib++)
   {
      const int SibPID = amr->patch[0][FaLv][FaPID]->sibling[sib];

      for (int d=0; d<3; d++)
      {
         loop      [d] = TABLE_01( sib, 'x'+d, FaGhost_Mag, PS1, FaGhost_Mag );
         offset_out[d] = TABLE_01( sib, 'x'+d, 0, FaGhost_Mag, FaGhost_Mag+PS1 );
      }

//    2.3.1 if the target sibling patch exists --> just copy data from it directly
      if ( SibPID >= 0 )
      {
         for (int d=0; d<3; d++)    offset_in[d] = TABLE_01( sib, 'x'+d, PS1-FaGhost_Mag, 0, 0 );

//       Bx
         if ( sib != 0  &&  sib != 1 ) // skip the normal direction
         {
            for (int k=0; k<loop[2]; k++)  {  k_out = k + offset_out[2];  k_in = k + offset_in[2];
            for (int j=0; j<loop[1]; j++)  {  j_out = j + offset_out[1];  j_in = j + offset_in[1];
                                              idx_B_in  = IDX321( 0, j_in,  k_in,  PS1P1,        PS1          );
                                              idx_B_out = IDX321( 0, j_out, k_out, FaSize_Mag_N, FaSize_Mag_T );
            for (int i=0; i<PS1P1;   i++)  {

               FaData_MagX[ idx_B_out ++ ] = amr->patch[FaSg_Mag][FaLv][SibPID]->magnetic[MAGX][ idx_B_in ++ ];

            }}}
         }

//       By
         if ( sib != 2  &&  sib != 3 ) // skip the normal direction
         {
            for (int k=0; k<loop[2]; k++)  {  k_out = k + offset_out[2];  k_in = k + offset_in[2];
            for (int j=0; j<PS1P1;   j++)  {  j_out = j;                  j_in = j;
                                              idx_B_in  = IDX321( offset_in[0],  j_in,  k_in,  PS1,          PS1P1        );
                                              idx_B_out = IDX321( offset_out[0], j_out, k_out, FaSize_Mag_T, FaSize_Mag_N );
            for (int i=0; i<loop[0]; i++)  {

               FaData_MagY[ idx_B_out ++ ] = amr->patch[FaSg_Mag][FaLv][SibPID]->magnetic[MAGY][ idx_B_in ++ ];

            }}}
         }

//       Bz
         if ( sib != 4  &&  sib != 5 ) // skip the normal direction
         {
            for (int k=0; k<PS1P1;   k++)  {  k_out = k;                  k_in = k;
            for (int j=0; j<loop[1]; j++)  {  j_out = j + offset_out[1];  j_in = j + offset_in[1];
                                              idx_B_in  = IDX321( offset_in[0],  j_in,  k_in,  PS1,          PS1          );
                                              idx_B_out = IDX321( offset_out[0], j_out, k_out, FaSize_Mag_T, FaSize_Mag_T );
            for (int i=0; i<loop[0]; i++)  {

               FaData_MagZ[ idx_B_out ++ ] = amr->patch[FaSg_Mag][FaLv][SibPID]->magnetic[MAGZ][ idx_B_in ++ ];

            }}}
         }
      } // if ( SibPID >= 0 )


//    2.3.2 if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
      else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
      {
//       work on one component at a time since the array sizes of different components are different
         for (int v=0; v<NCOMP_MAG; v++)
         {
//          get the normal direction
            const int norm_dir = ( v == MAGX ) ? 0 :
                                 ( v == MAGY ) ? 1 :
                                 ( v == MAGZ ) ? 2 : -1;
#           ifdef GAMER_DEBUG
            if ( norm_dir == -1 )   Aux_Error( ERROR_INFO, "Target face-centered variable != MAGX/Y/Z !!\n" );
#           endif

//          only need ghost zones along the two transverse directions
            if ( sib == norm_dir*2  ||  sib == norm_dir*2+1 )  continue;

//          set array indices --> correspond to the **cell-centered** array
            int FC_BC_Idx_Start[3], FC_BC_Idx_End[3], FC_BC_Size[3];
            double xyz_mag[3];   // cell-centered corner coordinates for the user-specified magnetic field B.C.
            for (int d=0; d<3; d++)
            {
               if ( d == norm_dir )
               {
                  FC_BC_Idx_Start[d] = 0;
                  FC_BC_Idx_End  [d] = FaSize_Mag_N - 2;
                  FC_BC_Size     [d] = FaSize_Mag_N - 1;
                  xyz_mag        [d] = amr->patch[0][FaLv][FaPID]->EdgeL[d] + 0.5*amr->dh[FaLv];
               }

               else
               {
                  FC_BC_Idx_Start[d] = offset_out[d];
                  FC_BC_Idx_End  [d] = loop[d] + FC_BC_Idx_Start[d] - 1;
                  FC_BC_Size     [d] = FaSize_Mag_T;
                  xyz_mag        [d] = amr->patch[0][FaLv][FaPID]->EdgeL[d] + (0.5-FaGhost_Mag)*amr->dh[FaLv];
               }
            }

            BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;
            real *FaData_Mag3v[NCOMP_MAG] = { FaData_MagX, FaData_MagY, FaData_MagZ };

            switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
            {
               case BC_FLU_OUTFLOW:
                  MHD_BoundaryCondition_Outflow   ( FaData_Mag3v, BC_Face[BC_Sibling], 1, FaGhost_Mag,
                                                    FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                    &v );
               break;

               case BC_FLU_REFLECTING:
                  MHD_BoundaryCondition_Reflecting( FaData_Mag3v, BC_Face[BC_Sibling], 1, FaGhost_Mag,
                                                    FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                    &v );
               break;

               case BC_FLU_DIODE:
                  MHD_BoundaryCondition_Diode     ( FaData_Mag3v, BC_Face[BC_Sibling], 1, FaGhost_Mag,
                                                    FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                    &v );
               break;

               case BC_FLU_USER:
                  MHD_BoundaryCondition_User      ( FaData_Mag3v, BC_Face[BC_Sibling], 1,
                                                    FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                    &v, Time[FaLv], amr->dh[FaLv], xyz_mag, FaLv );
               break;

               default:
                  Aux_Error( ERROR_INFO, "unsupported MHD B.C. (%d) !!\n", OPT__BC_FLU[ BC_Face[BC_Sibling] ] );

            } // switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
         } // for (int v=0; v<NCOMP_MAG; v++)
      } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )


//    2.3.3 it will violate the proper-nesting condition if the flagged patch is NOT surrounded by siblings
      else
         Aux_Error( ERROR_INFO, "SibPID = %d (FaLv %d, FaPID %d, Sib %d) !!\n", SibPID, FaLv, FaPID, sib );
   } // for (int sib=0; sib<6; sib++)
#  endif // #ifdef MHD

} // FUNCTION : PrepareCData



#endif // #ifdef LOAD_BALANCE
