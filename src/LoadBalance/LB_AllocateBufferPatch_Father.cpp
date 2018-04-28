#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_AllocateBufferPatch_Father
// Description :  Allocate father-buffer patches at SonLv-1 for the "real" patches at SonLv to ensure that
//                the proper-nesting condition is satisfied for these "real" patches
//
// Note        :  1. Eight father-buffer patches are allocated at a time
//                2. Father-buffer patches at SonLv-1 are NOT allocated for the "sibling-buffer" patches at SonLv
//                   --> Proper-nesting condition may not be satisfied for the "sibling-buffer" patches at SonLv
//                   --> However, all sibling-buffer patches at SonLv will have fathers after calling this function
//                       due to the proper-nesting condition of the real patches
//                       --> But note that, in the current implementation, the father indices of all sibling/father-buffer
//                           patches are always set to -1
//                3. Father-buffer patches at SonLv-1 are NOT allocated for the "father-buffer" patches at SonLv
//                4. This function will reconstruct the list LB_PaddedCr1DList[SonLv-1][]
//                5. SearchAllSon == true  --> search over all real patches at SonLv
//                                == false --> search over patches recorded in TargetSonPID0
//                6. RecordFaPID  == ture  --> record the indices of all newly-allocated father-buffer patches
//                                             (with LocalID==0)
//                                         --> the pointer "*NewFaBufPID0" will be allocated with memory
//                                         --> **this array must be deallocated manually**
//                                == false --> do nothing
//
// Parameter   :  SonLv         : Target refinement level of sons
//                SearchAllSon  : Whether to search over all real patches at SonLv or not
//                NInput        : Number of target son patches (with LocalID==0) in "TargetSonPID0"
//                                (useful only if "SearchAllSon == false")
//                TargetSonPID0 : Lists recording all target son patches (with LocalID==0)
//                                (useful only if "SearchAllSon == false")
//                RecordFaPID   : Record the indices of all newly-allocated father-buffer patches
//                                (with LocalID==0)
//                NNewFaBuf0    : Pointer recording the number of newly-allocated father-buffer patches
//                                (useful only if "RecordFaPID == true")
//                NewFaBufPID0  : Lists recording indices of all newly-allocated father-buffer patches
//                                with LocalID==0 (useful only if "RecordFaPID == true")
//-------------------------------------------------------------------------------------------------------
void LB_AllocateBufferPatch_Father( const int SonLv, const bool SearchAllSon, const int NInput, int* TargetSonPID0,
                                    const bool RecordFaPID, int* NNewFaBuf0, int** NewFaBufPID0 )
{

// base-level patches do not have fathers
   if ( SonLv == 0 )
   {
      if ( RecordFaPID )
      {
         *NNewFaBuf0   = 0;
         *NewFaBufPID0 = NULL;
      }

      return;
   }


// check
   if ( !SearchAllSon  &&  NInput != 0  &&  TargetSonPID0 == NULL )
      Aux_Error( ERROR_INFO, "SonLv %d, NInput %d, TargetSonPID0 == NULL !!\n", SonLv, NInput );

   if ( !SearchAllSon  &&  ( NInput < 0 || NInput > amr->NPatchComma[SonLv][1]/8 )  )
      Aux_Error( ERROR_INFO, "SonLv %d, incorrect NInput = %d (SonNReal0 = %d) !!\n",
                 SonLv, NInput, amr->NPatchComma[SonLv][1]/8 );

   if ( amr->NPatchComma[SonLv-1][3] != amr->num[SonLv-1] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][3] (%d) != amr->num[%d] (%d)\" !!\n",
                 SonLv-1, amr->NPatchComma[SonLv-1][3], SonLv-1, amr->num[SonLv-1] );


   const int   NTargetSon0         = ( SearchAllSon ) ? amr->NPatchComma[SonLv][1]/8 : NInput;
   const int   NFaBuf_Max          = NTargetSon0*8;    // check eight father patches for each son patch
   const int   FaLv                = SonLv-1;
   const int   Padded              = 1<<NLEVEL;
   const int   BoxNScale_Padded[3] = { amr->BoxScale[0]/PATCH_SIZE + 2*Padded,
                                       amr->BoxScale[1]/PATCH_SIZE + 2*Padded,
                                       amr->BoxScale[2]/PATCH_SIZE + 2*Padded }; //normalized and padded BoxScale
   const int   FaScale             = amr->scale[FaLv];
   const int   FaPScale            = PATCH_SIZE*FaScale;
   const int   FaPGScale           = 2*FaPScale;
   const ulong dr[3]               = { (ulong)FaScale,
                                       (ulong)FaScale*BoxNScale_Padded[0],
                                       (ulong)FaScale*BoxNScale_Padded[0]*BoxNScale_Padded[1] };
   const ulong dr2[3]              = { 2*dr[0], 2*dr[1], 2*dr[2] };
   const int   NP_Old              = amr->NPatchComma[FaLv][3];

// NFaBuf_Dup: # of father-buffer patches including the duplicated ones
   int   FaCr3D0[3], Start[3], FaCr3D[3], NFaBuf, NFaBuf_Dup, SonPID0;
   ulong FaCr1D0;

   ulong *FaCr1D_List = new ulong [NFaBuf_Max];


// nothing to do if there is no target real patch at SonLv
   if ( NTargetSon0 == 0 )
   {
      if ( RecordFaPID )
      {
         *NNewFaBuf0   = 0;
         *NewFaBufPID0 = NULL;
      }

      delete [] FaCr1D_List;

      return;
   }


// 0. construct the target son patch list
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


// 1. construct the candidate father-buffer patch list
   NFaBuf_Dup = 0;

   for (int t=0; t<NTargetSon0; t++)
   {
      SonPID0 = TargetSonPID0[t];
      FaCr1D0 = amr->patch[0][SonLv][SonPID0]->PaddedCr1D;

      for (int d=0; d<3; d++)
      {
         FaCr3D0[d] = amr->patch[0][SonLv][SonPID0]->corner[d];

         if ( FaCr3D0[d] % FaPGScale != 0 )
         {
            Start   [d]  = 0;
            FaCr3D0 [d] -= FaPScale;
            FaCr1D0     -= dr[d];

#           ifdef GAMER_DEBUG
            if ( FaCr1D0 < 0 )   Aux_Error( ERROR_INFO, "FaCr1D0 (%lu) < 0 (dr[%d] = %lu) !!\n", FaCr1D0, d, dr[d] );
#           endif
         }
         else
            Start[d] = -1;
      }

//    record the PaddedCr1D of the father-buffer patch candidates
      for (int k=Start[2]; k<=Start[2]+1; k++)
      {
//       no external buffer patches will be allocated for the non-periodic B.C.
         if ( OPT__BC_FLU[4] != BC_FLU_PERIODIC )
         {
            FaCr3D[2] = FaCr3D0[2] + k*FaPGScale;

            if ( FaCr3D[2] < 0  ||  FaCr3D[2] >= amr->BoxScale[2] )  continue;
         }

         for (int j=Start[1]; j<=Start[1]+1; j++)
         {
//          no external buffer patches will be allocated for the non-periodic B.C.
            if ( OPT__BC_FLU[2] != BC_FLU_PERIODIC )
            {
               FaCr3D[1] = FaCr3D0[1] + j*FaPGScale;

               if ( FaCr3D[1] < 0  ||  FaCr3D[1] >= amr->BoxScale[1] )  continue;
            }

            for (int i=Start[0]; i<=Start[0]+1; i++)
            {
//             no external buffer patches will be allocated for the non-periodic B.C.
               if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC )
               {
                  FaCr3D[0] = FaCr3D0[0] + i*FaPGScale;

                  if ( FaCr3D[0] < 0  ||  FaCr3D[0] >= amr->BoxScale[0] )  continue;
               }

#              ifdef GAMER_DEBUG
               if ( NFaBuf_Dup >= NFaBuf_Max )
                  Aux_Error( ERROR_INFO, "NFaBuf_Dup (%ld) >= NFaBuf_Max (%ld) !!\n", NFaBuf_Dup, NFaBuf_Max );
#              endif

               FaCr1D_List[ NFaBuf_Dup++ ] = FaCr1D0 + (ulong)i*dr2[0] + (ulong)j*dr2[1] + (ulong)k*dr2[2];

            } // i
         } // j
      } // k
   } // for (int t=0; t<NTargetSon0; t++)

// check
#  ifdef GAMER_DEBUG
   if (  OPT__BC_FLU[0] == BC_FLU_PERIODIC  &&  OPT__BC_FLU[2] == BC_FLU_PERIODIC  &&  OPT__BC_FLU[4] == BC_FLU_PERIODIC  &&
         NFaBuf_Dup != NFaBuf_Max  )
      Aux_Error( ERROR_INFO, "NFaBuf_Dup (%ld) != NFaBuf_Max (%ld) for the periodic BC !!\n", NFaBuf_Dup, NFaBuf_Max );
#  endif


// 2. sort the candidate list and remove duplicates (with the same FaCr1D)
   Mis_Heapsort( NFaBuf_Dup, FaCr1D_List, NULL );

   NFaBuf = ( NFaBuf_Dup > 0 ) ? 1 : 0;

   for (int t=1; t<NFaBuf_Dup; t++)
      if ( FaCr1D_List[t] != FaCr1D_List[t-1] )    FaCr1D_List[ NFaBuf++ ] = FaCr1D_List[t];


// 3. get the matching list
   char *Match = new char [NFaBuf];
   Mis_Matching_char( amr->num[FaLv], amr->LB->PaddedCr1DList[FaLv], NFaBuf, FaCr1D_List, Match );

#  ifdef GAMER_DEBUG
   if ( MPI_NRank == 1 )
   {
      for (int t=0; t<NFaBuf; t++)
      {
         if ( Match[t] != 1 )
            Aux_Error( ERROR_INFO, "FaCr1D_List[%d] = %lu has no matching father patch at level %d !!\n",
                       t, FaCr1D_List[t], FaLv );
      }
   }
#  endif


// 4. allocate the father-buffer patches
   int NNew0 = 0;

   if ( RecordFaPID )   *NewFaBufPID0 = new int [NFaBuf];

   for (int t=0; t<NFaBuf; t++)
   {
      if ( Match[t] != 1 )
      {
         Mis_Idx1D2Idx3D( BoxNScale_Padded, FaCr1D_List[t], FaCr3D );

         for (int d=0; d<3; d++)    FaCr3D[d] = ( FaCr3D[d] - Padded )*PATCH_SIZE;

         if ( RecordFaPID )   (*NewFaBufPID0)[ NNew0 ++ ] = amr->num[FaLv];

//       father patch is still unkown, data array is not allocated yet
         amr->pnew( FaLv, FaCr3D[0],          FaCr3D[1],          FaCr3D[2],          -1, false, false );
         amr->pnew( FaLv, FaCr3D[0]+FaPScale, FaCr3D[1],          FaCr3D[2],          -1, false, false );
         amr->pnew( FaLv, FaCr3D[0],          FaCr3D[1]+FaPScale, FaCr3D[2],          -1, false, false );
         amr->pnew( FaLv, FaCr3D[0],          FaCr3D[1],          FaCr3D[2]+FaPScale, -1, false, false );
         amr->pnew( FaLv, FaCr3D[0]+FaPScale, FaCr3D[1]+FaPScale, FaCr3D[2],          -1, false, false );
         amr->pnew( FaLv, FaCr3D[0],          FaCr3D[1]+FaPScale, FaCr3D[2]+FaPScale, -1, false, false );
         amr->pnew( FaLv, FaCr3D[0]+FaPScale, FaCr3D[1],          FaCr3D[2]+FaPScale, -1, false, false );
         amr->pnew( FaLv, FaCr3D[0]+FaPScale, FaCr3D[1]+FaPScale, FaCr3D[2]+FaPScale, -1, false, false );

//       check : no external buffer patches should be allocated for the non-periodic B.C.
#        ifdef GAMER_DEBUG
         for (int d=0; d<3; d++)
         {
            if (  OPT__BC_FLU[2*d] != BC_FLU_PERIODIC  &&  ( FaCr3D[d] < 0 || FaCr3D[d] >= amr->BoxScale[d] )  )
               Aux_Error( ERROR_INFO, "FaCr3D[%d] = %d lies outside the simulation box for non-periodic BC !!\n", d, FaCr3D[d] );
         }
#        endif

         amr->NPatchComma[FaLv][3] += 8;
      }
   } // for (int t=0; t<NFaBuf; t++)

   if ( RecordFaPID )   *NNewFaBuf0 = NNew0;

   for (int m=4; m<28; m++)   amr->NPatchComma[FaLv][m] = amr->NPatchComma[FaLv][3];

// check the amr->NPatchComma recording
   if ( amr->NPatchComma[FaLv][3] != amr->num[FaLv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][3] (%d) != amr->num[%d] (%d)\" !!\n",
                 FaLv, amr->NPatchComma[FaLv][3], FaLv, amr->num[FaLv] );


// 5. reconstruct LB_PaddedCr1DList and LB_PaddedCr1DList_Index_Table at SonLv-1
   const int NP_New = amr->NPatchComma[FaLv][3];

   if ( NP_New != NP_Old )
   {
      amr->LB->PaddedCr1DList         [FaLv] = (ulong*)realloc( amr->LB->PaddedCr1DList         [FaLv],
                                                                NP_New*sizeof(ulong) );
      amr->LB->PaddedCr1DList_IdxTable[FaLv] = (int*  )realloc( amr->LB->PaddedCr1DList_IdxTable[FaLv],
                                                                NP_New*sizeof(int ) );

//    must re-input all PaddedCr1D in order to get the correct LB_PaddedCr1DList_IdxTable[FaLv]
      for (int PID=0; PID<NP_New; PID++)
         amr->LB->PaddedCr1DList[FaLv][PID] = amr->patch[0][FaLv][PID]->PaddedCr1D;

      Mis_Heapsort( NP_New, amr->LB->PaddedCr1DList[FaLv], amr->LB->PaddedCr1DList_IdxTable[FaLv] );

//    check : no duplicate patches at FaLv
#     ifdef GAMER_DEBUG
      for (int t=1; t<amr->num[FaLv]; t++)
      {
         if ( amr->LB->PaddedCr1DList[FaLv][t] == amr->LB->PaddedCr1DList[FaLv][t-1] )
            Aux_Error( ERROR_INFO, "duplicate patches at lv %d, PaddedCr1D %lu, PID = %d and %d !!\n",
                       FaLv, amr->LB->PaddedCr1DList[FaLv][t], amr->LB->PaddedCr1DList_IdxTable[FaLv][t],
                       amr->LB->PaddedCr1DList_IdxTable[FaLv][t-1] );
      }
#     endif
   } // if ( NP_New != NP_Old )


// free memory
   delete [] FaCr1D_List;
   delete [] Match;
   if ( SearchAllSon )  delete [] TargetSonPID0;

} // FUNCTION : LB_AllocateBufferPatch_Father



#endif // #ifdef LOAD_BALANCE
