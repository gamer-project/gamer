#include "GAMER.h"

#ifdef LOAD_BALANCE



void PrepareCData( const int FaLv, const int FaPID, real *const FaData,
                   const int FaSg_Flu, const int FaGhost_Flu, const int NSide_Flu,
                   const int FaSg_Pot, const int FaGhost_Pot, const int NSide_Pot,
                   const int BC_Face[], const int FluVarIdxList[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Refine_GetNewRealPatchList
// Description :  Get the lists of father patches at FaLv to allocate/deallocate real son patches at FaLv+1
//
// Note        :  1. This function is invoked by the function "LB_Refine"
//                2. Coarse-grid data for creating new patches are also collected in "NewCData_Away"
//                   --> Data of all sibling-buffer patches at FaLv must be prepared in advance in order to
//                       prepare the coarse-grid data for spatial interpolation
//                3. Home/Away : target patches at home/not at home
//                4. Cr1D and CData lists are unsorted
//                5. Use "call-by-reference" for the input parameters
//
// Parameter   :  FaLv           : Target refinement level to be refined
//                NNew_Home      : Number of home patches at FaLv to allocate son patches
//                NewPID_Home    : Patch indices of home patches at FaLv to allocate son patches
//                NNew_Away      : Number of away patches at FaLv to allocate son patches
//                NewCr1D_Away   : Padded 1D corner of away patches at FaLv to allocate son patches
//                NewCData_Away  : Coarse-grid data of away patches at FaLv to allocate son patches
//                NDel_Home      : Number of home patches at FaLv to deallocate son patches
//                DelPID_Home    : Patch indices of home patches at FaLv to deallocate son patches
//                NDel_Away      : Number of away patches at FaLv to deallocate son patches
//                DelCr1D_Away   : Padded 1D corner of away patches at FaLv to deallocate son patches
//
//                PARTICLE-only parameters (call-by-reference)
//                RefineF2S_Send_NPatchTotal : Total number of patches for exchanging particles from fathers to sons
//                RefineF2S_Send_PIDList     : Patch indices for exchanging particles from fathers to sons
//                RefineF2S_Send_LBIdxList   : Load-balance indices for exchanging particles from fathers to sons
//
// Return      :  NNew_Home, NewPID_Home, NNew_Away, NewCr1D_Away, NewCData_Away, NDel_Home, DelPID_Home,
//                NDel_Away, DelCr1D_Away
//-------------------------------------------------------------------------------------------------------
void LB_Refine_GetNewRealPatchList( const int FaLv, int &NNew_Home, int *&NewPID_Home, int &NNew_Away,
                                    ulong *&NewCr1D_Away, real *&NewCData_Away, int &NDel_Home, int *&DelPID_Home,
                                    int &NDel_Away, ulong *&DelCr1D_Away,
                                    int &RefineF2S_Send_NPatchTotal, int *&RefineF2S_Send_PIDList,
                                    long *&RefineF2S_Send_LBIdxList )
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


// initialize variables
   for (int r=0; r<MPI_NRank; r++)
   {
      NewMemSize  [r] = MemUnit;
      DelMemSize  [r] = MemUnit;
      NewCr1D_Send[r] = (ulong*)malloc( NewMemSize[r]*sizeof(ulong) );
      DelCr1D_Send[r] = (ulong*)malloc( DelMemSize[r]*sizeof(ulong) );
      NewPID_Send [r] = (int*  )malloc( NewMemSize[r]*sizeof(int  ) );
      NNew_Send   [r] = 0;
      NDel_Send   [r] = 0;
   }
   NewPID_Home = (int*)malloc( NewMemSize[MPI_Rank]*sizeof(int) );
   DelPID_Home = (int*)malloc( DelMemSize[MPI_Rank]*sizeof(int) );
   NNew_Home   = 0;
   NDel_Home   = 0;


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


//       record the new lists
         if ( TRank == MPI_Rank ) // target son patches are home
         {
//          allocate enough memory
            if ( NNew_Home >= NewMemSize[TRank] )
            {
               NewMemSize[TRank] += MemUnit;
               NewPID_Home        = (int*)realloc( NewPID_Home, NewMemSize[TRank]*sizeof(int) );
            }

            NewPID_Home[ NNew_Home ++ ] = FaPID;
         }

         else // TRank != MPI_Rank (target son patches are not home)
         {
//          allocate enough memory
            if ( NNew_Send[TRank] >= NewMemSize[TRank] )
            {
               NewMemSize  [TRank] += MemUnit;
               NewCr1D_Send[TRank]  = (ulong*)realloc( NewCr1D_Send[TRank], NewMemSize[TRank]*sizeof(ulong) );
               NewPID_Send [TRank]  = (int*  )realloc( NewPID_Send [TRank], NewMemSize[TRank]*sizeof(int  ) );
            }

            NewCr1D_Send[TRank][ NNew_Send[TRank] ] = TP->PaddedCr1D;
            NewPID_Send [TRank][ NNew_Send[TRank] ] = FaPID;
            NNew_Send   [TRank] ++;

#           ifdef PARTICLE
#           ifdef DEBUG_PARTICLE
            if ( RefineF2S_Send_NPatchTotal >= FaNReal )
               Aux_Error( ERROR_INFO, "target index (%d) >= FaNReal (%d) !!\n", RefineF2S_Send_NPatchTotal, FaNReal );
#           endif

            RefineF2S_Send_PIDList  [RefineF2S_Send_NPatchTotal] = FaPID;
            RefineF2S_Send_LBIdxList[RefineF2S_Send_NPatchTotal] = LBIdx;  // this is the LBIdx of one of the sons

            RefineF2S_Send_NPatchTotal ++;
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
   int NSide_Flu, FaGhost_Flu;

   Int_Table( OPT__REF_FLU_INT_SCHEME, NSide_Flu, FaGhost_Flu );

   const int FaSize_Flu = PATCH_SIZE + 2*FaGhost_Flu;

#  ifdef GRAVITY
   const int FaSg_Pot = amr->PotSg[FaLv];
   int NSide_Pot, FaGhost_Pot;

   Int_Table( OPT__REF_POT_INT_SCHEME, NSide_Pot, FaGhost_Pot );

   const int FaSize_Pot = PATCH_SIZE + 2*FaGhost_Pot;
   const int PSize      = NCOMP_TOTAL*FaSize_Flu*FaSize_Flu*FaSize_Flu + FaSize_Pot*FaSize_Pot*FaSize_Pot;
#  else
   const int PSize      = NCOMP_TOTAL*FaSize_Flu*FaSize_Flu*FaSize_Flu;
#  endif

   int New_Send_Disp[MPI_NRank], New_Recv_Disp[MPI_NRank], NNew_Recv[MPI_NRank], NNew_Send_Total, NNew_Recv_Total;
   int Del_Send_Disp[MPI_NRank], Del_Recv_Disp[MPI_NRank], NDel_Recv[MPI_NRank], NDel_Send_Total, NDel_Recv_Total;
   int New_Send_Disp_CData[MPI_NRank], New_Recv_Disp_CData[MPI_NRank];
   int NNew_Send_CData[MPI_NRank], NNew_Recv_CData[MPI_NRank];
   int Counter;
   ulong *New_SendBuf_Cr1D=NULL, *New_RecvBuf_Cr1D=NULL, *Del_SendBuf_Cr1D=NULL, *Del_RecvBuf_Cr1D=NULL;
   real  *New_SendBuf_CData=NULL, *New_RecvBuf_CData=NULL;

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
      NNew_Send_CData    [r] = PSize*NNew_Send    [r];
      NNew_Recv_CData    [r] = PSize*NNew_Recv    [r];
      New_Send_Disp_CData[r] = PSize*New_Send_Disp[r];
      New_Recv_Disp_CData[r] = PSize*New_Recv_Disp[r];
   }

// variables to be returned
   NNew_Away         = NNew_Recv_Total;
   NDel_Away         = NDel_Recv_Total;
   NewCr1D_Away      = new ulong [NNew_Recv_Total      ];
   NewCData_Away     = new  real [NNew_Recv_Total*PSize];
   DelCr1D_Away      = new ulong [NDel_Recv_Total      ];

   New_SendBuf_Cr1D  = new ulong [NNew_Send_Total      ];
   New_RecvBuf_Cr1D  = NewCr1D_Away;
   New_SendBuf_CData = new real  [NNew_Send_Total*PSize];
   New_RecvBuf_CData = NewCData_Away;
   Del_SendBuf_Cr1D  = new ulong [NDel_Send_Total      ];
   Del_RecvBuf_Cr1D  = DelCr1D_Away;


// 2.3 prepare the MPI send buffers
// 2.3.1&2 new Cr1D/CData

// determine the priority of different boundary faces (z>y>x) to set the corner cells properly for the non-periodic B.C.
   int BC_Face[26], BC_Face_tmp[3], FluVarIdxList[NCOMP_TOTAL];

   for (int s=0; s<26; s++)
   {
      BC_Face_tmp[0] = TABLE_01( s, 'x', 0, -1, 1 );
      BC_Face_tmp[1] = TABLE_01( s, 'y', 2, -1, 3 );
      BC_Face_tmp[2] = TABLE_01( s, 'z', 4, -1, 5 );

//    z > y > x
      if      ( BC_Face_tmp[2] != -1 )   BC_Face[s] = BC_Face_tmp[2];
      else if ( BC_Face_tmp[1] != -1 )   BC_Face[s] = BC_Face_tmp[1];
      else if ( BC_Face_tmp[0] != -1 )   BC_Face[s] = BC_Face_tmp[0];
   }

   for (int v=0; v<NCOMP_TOTAL; v++)   FluVarIdxList[v] = v;

// prepare the coarse-grid data
   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<NNew_Send[r]; t++)
   {
      New_SendBuf_Cr1D[Counter] = NewCr1D_Send[r][t];

#     ifdef GRAVITY
      PrepareCData( FaLv, NewPID_Send[r][t], New_SendBuf_CData+Counter*PSize,
                    FaSg_Flu, FaGhost_Flu, NSide_Flu, FaSg_Pot, FaGhost_Pot, NSide_Pot, BC_Face, FluVarIdxList );
#     else
      PrepareCData( FaLv, NewPID_Send[r][t], New_SendBuf_CData+Counter*PSize,
                    FaSg_Flu, FaGhost_Flu, NSide_Flu, NULL_INT, NULL_INT, NULL_INT, BC_Face, FluVarIdxList );
#     endif

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

// 2.4.2 new CData
#  ifdef FLOAT8
   MPI_Alltoallv( New_SendBuf_CData, NNew_Send_CData, New_Send_Disp_CData, MPI_DOUBLE,
                  New_RecvBuf_CData, NNew_Recv_CData, New_Recv_Disp_CData, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( New_SendBuf_CData, NNew_Send_CData, New_Send_Disp_CData, MPI_FLOAT,
                  New_RecvBuf_CData, NNew_Recv_CData, New_Recv_Disp_CData, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif

// 2.4.3 delete Cr1D
   MPI_Alltoallv( Del_SendBuf_Cr1D, NDel_Send, Del_Send_Disp, MPI_UNSIGNED_LONG,
                  Del_RecvBuf_Cr1D, NDel_Recv, Del_Recv_Disp, MPI_UNSIGNED_LONG, MPI_COMM_WORLD );



// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( NewCr1D_Send[r] );
      free( DelCr1D_Send[r] );
      free( NewPID_Send [r] );
   }
   delete [] New_SendBuf_Cr1D;
   delete [] New_SendBuf_CData;
   delete [] Del_SendBuf_Cr1D;

} // FUNCTION : LB_Refine_GetNewRealPatchList



//-------------------------------------------------------------------------------------------------------
// Function    :  PrepareCData
// Description :  Prepare coarse-grid data for spatial interpolation
//
// Note        :  1. Data of all sibling-buffer patches at FaLv must be prepared in advance
//                2. This function is also used in "LB_Refine_AllocateNewPatch"
//
// Parameter   :  FaLv           : Coarse-grid refinement level
//                FaPID          : Father patch index to prepare the coarse-grid data
//                FaData         : Array to store the coarse-grid data
//                FaSg_Flu       : Sandglass for the fluid solver
//                FaGhost_Flu    : Ghost size for the fluid solver
//                NSide_Flu      : Number of sibling directions to prepare the ghost-zone data (6/26) for the fluid solver
//                FaSg_Pot       : Sandglass for the Poisson solver
//                FaGhost_Pot    : Ghost size for the Poisson solver
//                NSide_Pot      : Number of sibling directions to prepare the ghost-zone data (6/26) for the Poisson solver
//                BC_Face        : Corresponding boundary faces (0~5) along 26 sibling directions ->for non-periodic B.C. only
//                FluVarIdxList  : List of target fluid variable indices                          ->for non-periodic B.C. only
//
// Return      :  FaData
//-------------------------------------------------------------------------------------------------------
void PrepareCData( const int FaLv, const int FaPID, real *const FaData,
                   const int FaSg_Flu, const int FaGhost_Flu, const int NSide_Flu,
                   const int FaSg_Pot, const int FaGhost_Pot, const int NSide_Pot,
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
#  endif // #ifdef GAMER_DEBUG


// 1. fill up the central region of FaData
   const int FaSize_Flu   = PATCH_SIZE + 2*FaGhost_Flu;
   real *const FaData_Flu = FaData;
#  ifdef GRAVITY
   const int FaSize_Pot   = PATCH_SIZE + 2*FaGhost_Pot;
   real *const FaData_Pot = FaData + NCOMP_TOTAL*FaSize_Flu*FaSize_Flu*FaSize_Flu;
#  endif
   int Idx, I, J, K;

// 1.1 fluid data
   for (int v=0; v<NCOMP_TOTAL; v++)   {
   for (int k=0; k<PATCH_SIZE; k++)    {  K = k + FaGhost_Flu;
   for (int j=0; j<PATCH_SIZE; j++)    {  J = j + FaGhost_Flu;
   for (int i=0; i<PATCH_SIZE; i++)    {  I = i + FaGhost_Flu;

      Idx = ((v*FaSize_Flu + K)*FaSize_Flu + J)*FaSize_Flu + I;

      FaData_Flu[Idx] = amr->patch[FaSg_Flu][FaLv][FaPID]->fluid[v][k][j][i];

   }}}}

// 1.2 potential data
#  ifdef GRAVITY
   for (int k=0; k<PATCH_SIZE; k++)    {  K = k + FaGhost_Pot;
   for (int j=0; j<PATCH_SIZE; j++)    {  J = j + FaGhost_Pot;
   for (int i=0; i<PATCH_SIZE; i++)    {  I = i + FaGhost_Pot;

      Idx = (K*FaSize_Pot + J)*FaSize_Pot + I;

      FaData_Pot[Idx] = amr->patch[FaSg_Pot][FaLv][FaPID]->pot[k][j][i];

   }}}
#  endif


// 2. fill up the ghost zone of FaData (no interpolation is required)
   const int  NDer       = 0;
   const int *DerVarList = NULL;

   int    Loop[3], Disp1[3], Disp2[3], I2, J2, K2, SibPID;
   int    BC_Sibling, BC_Idx_Start[3], BC_Idx_End[3];
   double xyz[3];

// calculate the corner coordinates of the coarse-grid data for the user-specified B.C.
   for (int d=0; d<3; d++)    xyz[d] = amr->patch[0][FaLv][FaPID]->EdgeL[d] + (0.5-FaGhost_Flu)*amr->dh[FaLv];

// 2.1 fluid data
   for (int sib=0; sib<NSide_Flu; sib++)
   {
      SibPID = amr->patch[0][FaLv][FaPID]->sibling[sib];

      for (int d=0; d<3; d++)
      {
         Loop [d] = TABLE_01( sib, 'x'+d, FaGhost_Flu, PATCH_SIZE, FaGhost_Flu );
         Disp1[d] = TABLE_01( sib, 'x'+d, 0, FaGhost_Flu, FaGhost_Flu+PATCH_SIZE );
      }

//    2.1.1 if the target sibling patch exists --> just copy data from the nearby patches at the same level
      if ( SibPID >= 0 )
      {
         for (int d=0; d<3; d++)    Disp2[d] = TABLE_01( sib, 'x'+d, PATCH_SIZE-FaGhost_Flu, 0, 0 );

         for (int v=0; v<NCOMP_TOTAL; v++){
         for (int k=0; k<Loop[2]; k++)    {  K = k + Disp1[2];    K2 = k + Disp2[2];
         for (int j=0; j<Loop[1]; j++)    {  J = j + Disp1[1];    J2 = j + Disp2[1];
         for (int i=0; i<Loop[0]; i++)    {  I = i + Disp1[0];    I2 = i + Disp2[0];

            Idx = ((v*FaSize_Flu + K)*FaSize_Flu + J)*FaSize_Flu + I;

            FaData_Flu[Idx] = amr->patch[FaSg_Flu][FaLv][SibPID]->fluid[v][K2][J2][I2];

         }}}}
      }


//    2.1.2 if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
      else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
      {
         for (int d=0; d<3; d++)
         {
            BC_Idx_Start[d] = Disp1[d];
            BC_Idx_End  [d] = Loop[d] + BC_Idx_Start[d] - 1;
         }

         BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;

         switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
         {
            case BC_FLU_OUTFLOW:
               Hydro_BoundaryCondition_Outflow   ( FaData_Flu, BC_Face[BC_Sibling], NCOMP_TOTAL, FaGhost_Flu,
                                                   FaSize_Flu, FaSize_Flu, FaSize_Flu, BC_Idx_Start, BC_Idx_End );
            break;

#           if ( MODEL == HYDRO  ||  MODEL == MHD )
            case BC_FLU_REFLECTING:
               Hydro_BoundaryCondition_Reflecting( FaData_Flu, BC_Face[BC_Sibling], NCOMP_TOTAL, FaGhost_Flu,
                                                   FaSize_Flu, FaSize_Flu, FaSize_Flu, BC_Idx_Start, BC_Idx_End,
                                                   FluVarIdxList, NDer, DerVarList );
            break;
#           if ( MODEL == MHD )
#           warning : WAIT MHD !!!
#           endif
#           endif

            case BC_FLU_USER:
               Flu_BoundaryCondition_User        ( FaData_Flu,                      NCOMP_TOTAL,
                                                   FaSize_Flu, FaSize_Flu, FaSize_Flu, BC_Idx_Start, BC_Idx_End,
                                                   FluVarIdxList, Time[FaLv], amr->dh[FaLv], xyz, _TOTAL, FaLv );
            break;

            default:
               Aux_Error( ERROR_INFO, "unsupported fluid B.C. (%d) !!\n", OPT__BC_FLU[ BC_Face[BC_Sibling] ] );

         } // switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
      } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )


//    2.1.3 it will violate the proper-nesting condition if the flagged patch is NOT surrounded by siblings
      else if ( SibPID == -1 )
         Aux_Error( ERROR_INFO, "no sibling patch is found for FaLv %d, FaPID %d, sib %d !!\n", FaLv, FaPID, sib );

      else
         Aux_Error( ERROR_INFO, "SibPID == %d (FaPID %d, sib %d) !!\n", SibPID, FaPID, sib );

   } // for (int sib=0; sib<NSide_Flu; sib++)


// 2.2 potential data
#  ifdef GRAVITY
   for (int sib=0; sib<NSide_Pot; sib++)
   {
      SibPID = amr->patch[0][FaLv][FaPID]->sibling[sib];

//    it will violate the proper-nesting condition if the flagged patch is NOT surrounded by sibling patches
#     ifdef GAMER_DEBUG
      if ( SibPID < 0 )
         Aux_Error( ERROR_INFO, "no sibling patch is found for FaLv %d, FaPID %d, sib %d !!\n", FaLv, FaPID, sib );
#     endif

      for (int d=0; d<3; d++)
      {
         Loop [d] = TABLE_01( sib, 'x'+d, FaGhost_Pot, PATCH_SIZE, FaGhost_Pot );
         Disp1[d] = TABLE_01( sib, 'x'+d, 0, FaGhost_Pot, FaGhost_Pot+PATCH_SIZE );
         Disp2[d] = TABLE_01( sib, 'x'+d, PATCH_SIZE-FaGhost_Pot, 0, 0 );
      }

      for (int k=0; k<Loop[2]; k++)    {  K = k + Disp1[2];    K2 = k + Disp2[2];
      for (int j=0; j<Loop[1]; j++)    {  J = j + Disp1[1];    J2 = j + Disp2[1];
      for (int i=0; i<Loop[0]; i++)    {  I = i + Disp1[0];    I2 = i + Disp2[0];

         Idx = (K*FaSize_Pot + J)*FaSize_Pot + I;

         FaData_Pot[Idx] = amr->patch[FaSg_Pot][FaLv][SibPID]->pot[K2][J2][I2];
      }}}

   } // for (int sib=0; sib<NSide_Pot; sib++)
#  endif // #ifdef GRAVITY

} // FUNCTION : PrepareCData



#endif // #ifdef LOAD_BALANCE
